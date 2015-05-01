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

#include <stk_mesh/base/MetaData.hpp>
#include <string.h>                     // for strcmp, strncmp
#include <Shards_CellTopologyManagedData.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <set>                          // for set
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CellTopology.hpp>  // for CellTopology
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceMin, etc
#include <stk_util/util/string_case_compare.hpp>  // for equal_case
#include "Shards_BasicTopologies.hpp"   // for getCellTopologyData, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part, etc
#include "stk_mesh/base/PropertyBase.hpp"  // for Property
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_mesh/baseImpl/PartRepository.hpp"  // for PartRepository
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.tcc"    // for topology::num_nodes, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc

namespace stk {
namespace mesh {

namespace {

bool root_part_in_subset(stk::mesh::Part & part)
{
  if (is_topology_root_part(part)) {
    return true;
  }
  const PartVector & subsets = part.subsets();
  for (PartVector::const_iterator it=subsets.begin() ; it != subsets.end() ; ++it) {
    if (is_topology_root_part( **it )) {
      return true;
    }
  }
  return false;
}

void find_cell_topologies_in_part_and_subsets_of_same_rank(const Part & part, EntityRank rank, std::set<CellTopology> & topologies_found)
{
  MetaData & meta = MetaData::get(part);
  CellTopology top = meta.get_cell_topology(part);
  if ((top.isValid() && (part.primary_entity_rank() == rank))) {
    topologies_found.insert(top);
  }
  const PartVector & subsets = part.subsets();
  for (PartVector::const_iterator it=subsets.begin() ; it != subsets.end() ; ++it) {
    top = meta.get_cell_topology(**it);
    if (top.isValid() && ( (**it).primary_entity_rank() == rank) ) {
      topologies_found.insert(top);
    }
  }
}

//----------------------------------------------------------------------

stk::mesh::FieldBase* try_to_find_coord_field(const stk::mesh::MetaData& meta)
{
  //attempt to initialize the coordinate-field pointer, trying a couple
  //of commonly-used names. It is expected that the client code will initialize
  //the coordinates field using set_coordinate_field, but this is an
  //attempt to be helpful for existing client codes which aren't yet calling that.

  stk::mesh::FieldBase* coord_field = meta.get_field(stk::topology::NODE_RANK, "mesh_model_coordinates");
  if (coord_field == NULL) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "mesh_model_coordinates_0");
  }
  if (coord_field == NULL) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "model_coordinates");
  }
  if (coord_field == NULL) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "model_coordinates_0");
  }
  if (coord_field == NULL) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "coordinates");
  }

  return coord_field;
}

} // namespace

void MetaData::assign_cell_topology(
  Part &                   part,
  const CellTopology       cell_topology)
{
  const size_t part_ordinal = part.mesh_meta_data_ordinal();

  if (part_ordinal >= m_partCellTopologyVector.size()) {
    m_partCellTopologyVector.resize(part_ordinal + 1);
  }

  m_partCellTopologyVector[part_ordinal] = cell_topology;

  stk::topology topo = stk::mesh::get_topology(cell_topology, m_spatial_dimension);

  m_part_repo.get_all_parts()[part_ordinal]->m_partImpl.set_topology(topo);

  ThrowRequireMsg(cell_topology.getCellTopologyData(), "bad topology in MetaData::assign_cell_topology");
}

void MetaData::set_mesh_on_fields(BulkData* bulk)
{
  const FieldVector& fields = get_fields();
  for(size_t i=0; i<fields.size(); ++i) {
    fields[i]->set_mesh(bulk);
  }
}

MetaData & MetaData::get( const BulkData & bulk_data) {
  return bulk_data.meta_data();
}

MetaData & MetaData::get( const Bucket & bucket) {
  return MetaData::get(BulkData::get(bucket));
}

MetaData & MetaData::get( const Ghosting & ghost) {
  return MetaData::get(BulkData::get(ghost));
}
//----------------------------------------------------------------------

std::ostream &
print_entity_id( std::ostream & os , const MetaData & meta_data ,
                  EntityRank type , EntityId id )
{
  const std::string & name = meta_data.entity_rank_name( type );
  return os << name << "[" << id << "]" ;
}


std::ostream &
print_entity_key( std::ostream & os , const MetaData & meta_data ,
                  const EntityKey & key )
{
  const EntityRank type   = key.rank();
  const EntityId   id     = key.id();
  return print_entity_id( os , meta_data , type , id );
}

std::string
print_entity_key( const MetaData & meta_data , const EntityKey & key )
{
  std::ostringstream out;
  print_entity_key(out, meta_data, key);
  return out.str();
}

//----------------------------------------------------------------------

void MetaData::require_not_committed() const
{
  ThrowRequireMsg(!m_commit, "mesh MetaData has been committed.");
}

void MetaData::require_committed() const
{
  ThrowRequireMsg(m_commit, "mesh MetaData has not been committed.");
}

void MetaData::require_same_mesh_meta_data( const MetaData & rhs ) const
{
  ThrowRequireMsg(this == &rhs, "Different mesh_meta_data.");
}

void MetaData::require_valid_entity_rank( EntityRank rank ) const
{
  ThrowRequireMsg(check_rank(rank),
      "entity_rank " << rank << " >= " << m_entity_rank_names.size() );
  ThrowRequireMsg( !(rank == stk::topology::FACE_RANK && spatial_dimension() == 2),
                   "Should not use FACE_RANK in 2d");
}

void MetaData::require_not_relation_target( const Part * const part ) const
{

}

//----------------------------------------------------------------------

MetaData::MetaData(size_t spatial_dimension, const std::vector<std::string>& entity_rank_names)
  : m_bulk_data(NULL),
    m_commit( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_aura_part(NULL),
    m_field_repo(),
    m_coord_field(NULL),
    m_properties( ),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/),
    m_side_rank(stk::topology::INVALID_RANK),
    m_part_fields()
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_internal_part("OWNS");
  m_shares_part = & declare_internal_part("SHARES");
  m_aura_part = & declare_internal_part("AURA");

  initialize(spatial_dimension, entity_rank_names);
}

MetaData::MetaData()
  : m_bulk_data(NULL),
    m_commit( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_aura_part(NULL),
    m_field_repo(),
    m_coord_field(NULL),
    m_properties( ),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/),
    m_side_rank(stk::topology::INVALID_RANK),
    m_part_fields()
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_internal_part("OWNS");
  m_shares_part = & declare_internal_part("SHARES");
  m_aura_part = & declare_internal_part("AURA");
}

//----------------------------------------------------------------------

void MetaData::initialize(size_t spatial_dimension, const std::vector<std::string> &rank_names)
{
  ThrowErrorMsgIf( !m_entity_rank_names.empty(), "already initialized");
  ThrowErrorMsgIf( spatial_dimension > 3, "Max spatial dimension is 3");

  if ( rank_names.empty() ) {
    m_entity_rank_names = stk::mesh::entity_rank_names();
  }
  else {
    ThrowErrorMsgIf(rank_names.size() < stk::topology::ELEMENT_RANK+1,
                    "Entity rank name vector must name every rank, rank_names.size() = " <<
                    rank_names.size() << ", need " << stk::topology::ELEMENT_RANK+1 << " names");
    m_entity_rank_names = rank_names;
  }

  m_spatial_dimension = spatial_dimension;
  m_side_rank = side_rank();

  internal_declare_known_cell_topology_parts();
}

const std::string& MetaData::entity_rank_name( EntityRank entity_rank ) const
{
  ThrowErrorMsgIf( entity_rank >= m_entity_rank_names.size(),
      "entity-rank " << entity_rank <<
      " out of range. Must be in range 0.." << m_entity_rank_names.size());

  return m_entity_rank_names[entity_rank];
}

EntityRank MetaData::entity_rank( const std::string &name ) const
{
  EntityRank entity_rank = InvalidEntityRank;

  for (size_t i = 0; i < m_entity_rank_names.size(); ++i)
    if (equal_case(name, m_entity_rank_names[i])) {
      entity_rank = static_cast<EntityRank>(i);
      break;
    }
  return entity_rank;
}

FieldBase const* MetaData::coordinate_field() const
{
  if (m_coord_field == NULL) {
    m_coord_field = try_to_find_coord_field(*this);
  }

  ThrowErrorMsgIf( m_coord_field == NULL,
                   "MetaData::coordinate_field: Coordinate field has not been defined" );

  return m_coord_field;
}

//----------------------------------------------------------------------

Part * MetaData::get_part( const std::string & p_name ,
                           const char * required_by ) const
{
  const PartVector & all_parts = m_part_repo.get_all_parts();

  Part * const p = find( all_parts , p_name );

  ThrowErrorMsgIf( required_by && NULL == p,
                   "Failed to find part with name " << p_name <<
                   " for method " << required_by );

  return p ;
}

void MetaData::add_new_part_in_part_fields()
{
  for(size_t i=0; i<m_part_fields.size(); ++i) {
    PartFieldBase* part_field = m_part_fields[i];
    std::vector<char*>& char_ptr_vector = part_field->char_data();
    char_ptr_vector.push_back(new char[part_field->bytes_per_part()]);
  }
}

void MetaData::synchronize_part_fields_with_parts()
{
  size_t num_parts = get_parts().size();
  for(size_t i=0; i<m_part_fields.size(); ++i) {
    PartFieldBase* part_field = m_part_fields[i];
    std::vector<char*>& char_ptr_vector = part_field->char_data();
    if (char_ptr_vector.size() != num_parts) {
      size_t old_size = char_ptr_vector.size();
      char_ptr_vector.resize(num_parts);
      for(size_t j=old_size; j<num_parts; ++j) {
        char_ptr_vector[j] = new char[part_field->bytes_per_part()];
      }
    }
  }
}

Part & MetaData::declare_part( const std::string & p_name )
{
  const EntityRank rank = InvalidEntityRank;

  add_new_part_in_part_fields();
  return *m_part_repo.declare_part( p_name, rank );
}

Part & MetaData::declare_internal_part( const std::string & p_name )
{
  std::string internal_name = impl::convert_to_internal_name(p_name);
  return declare_part(internal_name);
}

Part & MetaData::declare_part( const std::string & p_name , EntityRank rank, bool arg_force_no_induce )
{
  require_valid_entity_rank(rank);

  add_new_part_in_part_fields();
  return *m_part_repo.declare_part( p_name , rank, arg_force_no_induce );
}

Part & MetaData::declare_internal_part( const std::string & p_name , EntityRank rank )
{
  std::string internal_name = impl::convert_to_internal_name(p_name);
  return declare_part(internal_name, rank);
}

void MetaData::declare_part_subset( Part & superset , Part & subset )
{
  if (!is_initialized()) {
    // can't do any topology stuff yet
    return internal_declare_part_subset(superset, subset);
  }

  CellTopology superset_top = get_cell_topology(superset);

  const bool no_superset_topology = !superset_top.isValid();
  if ( no_superset_topology ) {
    internal_declare_part_subset(superset,subset);
    return;
  }
  // Check for cell topology root parts in subset or subset's subsets
  const bool subset_has_root_part = root_part_in_subset(subset);
  ThrowErrorMsgIf( subset_has_root_part, "MetaData::declare_part_subset:  Error, root cell topology part found in subset or below." );

  std::set<CellTopology> cell_topologies;
  find_cell_topologies_in_part_and_subsets_of_same_rank(subset,superset.primary_entity_rank(),cell_topologies);

  ThrowErrorMsgIf( cell_topologies.size() > 1,
      "MetaData::declare_part_subset:  Error, multiple cell topologies of rank "
      << superset.primary_entity_rank()
      << " defined below subset"
      );
  const bool non_matching_cell_topology = ((cell_topologies.size() == 1) && (*cell_topologies.begin() != superset_top));
  ThrowErrorMsgIf( non_matching_cell_topology,
      "MetaData::declare_part_subset:  Error, superset topology = "
      << superset_top.getName() << " does not match the topology = "
      << cell_topologies.begin()->getName()
      << " coming from the subset part"
      );
  // Everything is Okay!
  internal_declare_part_subset(superset,subset);
  // Update PartCellTopologyVector for "subset" and same-rank subsets, ad nauseum
  if (subset.primary_entity_rank() == superset.primary_entity_rank()) {
    assign_cell_topology( subset, superset_top);
    const PartVector & subset_parts = subset.subsets();
    for (PartVector::const_iterator it=subset_parts.begin() ; it != subset_parts.end() ; ++it) {
      Part & it_part = **it;
      if (it_part.primary_entity_rank() == superset.primary_entity_rank()) {
        assign_cell_topology( it_part, superset_top);
      }
    }
  }
}

void MetaData::internal_declare_part_subset( Part & superset , Part & subset )
{
  require_not_committed();
  require_same_mesh_meta_data( MetaData::get(superset) );
  require_same_mesh_meta_data( MetaData::get(subset) );
  require_not_relation_target( &superset );
  require_not_relation_target( &subset );

  m_part_repo.declare_subset( superset, subset );

  // The new superset / subset relationship can cause a
  // field restriction to become incompatible or redundant.
  m_field_repo.verify_and_clean_restrictions(superset, subset);
}

//----------------------------------------------------------------------

void MetaData::declare_field_restriction(
  FieldBase      & arg_field ,
  const Part     & arg_part ,
  const unsigned   arg_num_scalars_per_entity ,
  const unsigned   arg_first_dimension ,
  const void     * arg_init_value )
{
  static const char method[] =
    "std::mesh::MetaData::declare_field_restriction" ;

  require_same_mesh_meta_data( MetaData::get(arg_field) );
  require_same_mesh_meta_data( MetaData::get(arg_part) );

  m_field_repo.declare_field_restriction(
      method,
      arg_field,
      arg_part,
      m_part_repo.get_all_parts(),
      arg_num_scalars_per_entity,
      arg_first_dimension,
      arg_init_value
      );
}

void MetaData::declare_field_restriction(
  FieldBase      & arg_field ,
  const Selector & arg_selector ,
  const unsigned   arg_num_scalars_per_entity ,
  const unsigned   arg_first_dimension ,
  const void     * arg_init_value )
{
  static const char method[] =
    "std::mesh::MetaData::declare_field_restriction" ;

  require_same_mesh_meta_data( MetaData::get(arg_field) );

  m_field_repo.declare_field_restriction(
      method,
      arg_field,
      arg_selector,
      m_part_repo.get_all_parts(),
      arg_num_scalars_per_entity,
      arg_first_dimension,
      arg_init_value
      );
}

//----------------------------------------------------------------------

void MetaData::commit()
{
  require_not_committed();

  m_commit = true ; // Cannot add or change parts or fields now

  synchronize_part_fields_with_parts();

  set_mesh_on_fields(m_bulk_data);

#ifdef STK_VERBOSE_OUTPUT
  dump_all_meta_info(std::cout);
#endif
}

MetaData::~MetaData()
{
  // Destroy part fields:
  for(size_t i=0; i<m_part_fields.size(); ++i) {
    delete m_part_fields[i];
  }

  // Destroy the properties, used 'new' to allocate so now use 'delete'

  try {
    std::vector<PropertyBase * >::iterator j = m_properties.begin();

    for ( ; j != m_properties.end() ; ++j ) { delete *j ; }

    m_properties.clear();

    std::vector<shards::CellTopologyManagedData*>::iterator i = m_created_topologies.begin();
    for ( ; i != m_created_topologies.end(); ++i) {
      delete *i;
    }
  } catch(...) {}

  // PartRepository is member data
  // FieldRepository is member data
}

void MetaData::internal_declare_known_cell_topology_parts()
{
  // Load up appropriate standard cell topologies.
  register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Node >()), stk::topology::NODE_RANK);

  if (m_spatial_dimension == 1) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), stk::topology::ELEMENT_RANK); // ???
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), stk::topology::ELEMENT_RANK); // ???

  }

  else if (m_spatial_dimension == 2) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()), stk::topology::ELEMENT_RANK);
  }

  else if (m_spatial_dimension == 3) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), stk::topology::EDGE_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), stk::topology::EDGE_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<11> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()), stk::topology::ELEMENT_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()), stk::topology::ELEMENT_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()), stk::topology::ELEMENT_RANK);
  }
}

void MetaData::register_cell_topology(const CellTopology cell_topology, EntityRank entity_rank)
{
  ThrowRequireMsg(is_initialized(),"MetaData::register_cell_topology: initialize() must be called before this function");

  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);

  const bool       duplicate     = it != m_cellTopologyPartEntityRankMap.end();
  const EntityRank existing_rank = duplicate ? (*it).second.second : stk::topology::NODE_RANK;

  ThrowErrorMsgIf(duplicate && existing_rank != entity_rank,
    "For args: cell_topolgy " << cell_topology.getName() << " and entity_rank " << entity_rank << ", " <<
    "previously declared rank = " << existing_rank );

  if (! duplicate) {
    std::string part_name = std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + std::string(cell_topology.getName());

    ThrowErrorMsgIf(get_part(part_name) != 0, "Cannot register topology with same name as existing part '" << cell_topology.getName() << "'" );

    Part &part = declare_internal_part(part_name, entity_rank);
    m_cellTopologyPartEntityRankMap[cell_topology] = CellTopologyPartEntityRankMap::mapped_type(&part, entity_rank);

    assign_cell_topology( part, cell_topology);
  }
  //check_topo_db();
}

shards::CellTopology MetaData::register_superelement_cell_topology(stk::topology topo)
{
  shards::CellTopology cell_topology = get_cell_topology(topo.name());
  if (!cell_topology.isValid()) {
    shards::CellTopologyManagedData *cell_topology_data = new shards::CellTopologyManagedData(topo.name());
    m_created_topologies.push_back(cell_topology_data);
    cell_topology = shards::CellTopology(cell_topology_data);

    cell_topology_data->base              = cell_topology_data ;
    cell_topology_data->dimension         = 1 ;
    cell_topology_data->vertex_count      = topo.num_nodes();
    cell_topology_data->node_count        = topo.num_nodes();
    cell_topology_data->edge_count        = 0 ;
    cell_topology_data->side_count        = 0 ;
    cell_topology_data->permutation_count = 0 ;
    cell_topology_data->subcell_count[0]  = topo.num_nodes();
    cell_topology_data->subcell_count[1]  = 1 ;
    cell_topology_data->subcell_count[2]  = 0 ;
    cell_topology_data->subcell_count[3]  = 0 ;

    register_cell_topology(cell_topology, stk::topology::ELEMENT_RANK);
  }
  return cell_topology;
}

CellTopology
MetaData::get_cell_topology(
  const std::string &   topology_name) const
{
  std::string part_name = impl::convert_to_internal_name(std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + topology_name);

  Part *part = get_part(part_name);
  if (part)
    return get_cell_topology(*part);
  else
    return CellTopology();
}

Part &MetaData::get_cell_topology_root_part(const CellTopology cell_topology) const
{
  ThrowRequireMsg(is_initialized(),"MetaData::get_cell_topology_root_part: initialize() must be called before this function");
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  ThrowErrorMsgIf(it == m_cellTopologyPartEntityRankMap.end(),
                  "Cell topology " << cell_topology.getName() <<
                  " has not been registered");

  return *(*it).second.first;
}

/// Note:  This function only uses the PartCellTopologyVector to look up the
/// cell topology for a given part.
/// This depends on declare_part_subset to update this vector correctly.  If a
/// cell topology is not defined for the given part, then an invalid Cell
/// Topology object will be returned.
CellTopology MetaData::get_cell_topology( const Part & part) const
{
  ThrowRequireMsg(is_initialized(),"MetaData::get_cell_topology: initialize() must be called before this function");
  CellTopology cell_topology;

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  if (part_ordinal < m_partCellTopologyVector.size())
    {
      cell_topology = m_partCellTopologyVector[part_ordinal];
    }

  return cell_topology;
}

stk::topology MetaData::get_topology(const Part & part) const
{
    return stk::mesh::get_topology(get_cell_topology(part), spatial_dimension());
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Verify parallel consistency of fields and parts

namespace {

void pack( CommBuffer & b , const PartVector & pset )
{
  PartVector::const_iterator i , j ;
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartVector & subsets   = p.subsets();

    const size_t       name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    {
      const unsigned ord = p.mesh_meta_data_ordinal();
      b.pack<unsigned>( ord );
    }

    b.pack<unsigned>( name_len );
    b.pack<char>( name_ptr , name_len );

    const unsigned subset_size = static_cast<unsigned>(subsets.size());
    b.pack<unsigned>( subset_size );
    for ( j = subsets.begin() ; j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      const unsigned ord = s.mesh_meta_data_ordinal();
      b.pack<unsigned>( ord );
    }
  }
}

bool unpack_verify( CommBuffer & b , const PartVector & pset )
{
  enum { MAX_TEXT_LEN = 4096 };
  char b_text[ MAX_TEXT_LEN ];
  unsigned b_tmp = 0;

  bool ok = true ;
  PartVector::const_iterator i , j ;
  for ( i = pset.begin() ; ok && i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartVector & subsets   = p.subsets();
    const unsigned     name_len = static_cast<unsigned>(p.name().size()) + 1 ;
    const char * const name_ptr = p.name().c_str();

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == p.mesh_meta_data_ordinal();
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == name_len ;
    }
    if ( ok ) {
      b.unpack<char>( b_text , name_len );
      ok = 0 == strcmp( name_ptr , b_text );
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == subsets.size() ;
    }
    for ( j = subsets.begin() ; ok && j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == s.mesh_meta_data_ordinal();
    }

  }
  return ok ;
}

void pack( CommBuffer & ,
           const std::vector< FieldBase * > & )
{
}

bool unpack_verify( CommBuffer & ,
                    const std::vector< FieldBase * > & )
{
  bool ok = true ;
  return ok ;
}

}

//----------------------------------------------------------------------

void verify_parallel_consistency( const MetaData & s , ParallelMachine pm )
{
  const unsigned p_rank = parallel_machine_rank( pm );

  const bool is_root = 0 == p_rank ;

  CommBroadcast comm( pm , 0 );

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields() );
  }

  comm.allocate_buffer();

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields() );
  }

  comm.communicate();

  int ok[ 2 ];

  ok[0] = unpack_verify( comm.recv_buffer() , s.get_parts() );
  ok[1] = unpack_verify( comm.recv_buffer() , s.get_fields() );

  all_reduce( pm , ReduceMin<2>( ok ) );

  ThrowRequireMsg(ok[0], "P" << p_rank << ": FAILED for Parts");
  ThrowRequireMsg(ok[1], "P" << p_rank << ": FAILED for Fields");
}

//----------------------------------------------------------------------

bool is_topology_root_part(const Part & part) {
  MetaData & meta = MetaData::get(part);
  CellTopology top = meta.get_cell_topology(part);
  if (top.isValid()) {
    const Part & root_part = meta.get_cell_topology_root_part(top);
    return (root_part == part);
  }
  return false;
}

/// This is a convenience function to get the root cell topology part and then
/// call declare_part_subset.
/// Note:  MetaData::declare_part_subset is the function that actually
/// updates the PartCellTopologyVector in MetaData for fast look-up of the
/// Cell Topology.
void set_topology(Part & part, stk::topology topo)
{
  MetaData& meta = part.mesh_meta_data();
  if (part.primary_entity_rank() == InvalidEntityRank) {
    //declare_part will set the rank on the part, if the part already exists.
    meta.declare_part(part.name(), topo.rank());
  }

  if (topo.is_superelement()) {
    // Need to (possibly) create a CellTopology corresponding to this superelement stk::topology.
    shards::CellTopology cell_topology = meta.register_superelement_cell_topology(topo);
    set_cell_topology(part, cell_topology);
  } else {
    set_cell_topology(part, get_cell_topology(topo));
  }
}

void set_cell_topology(
  Part &                        part,
  CellTopology             cell_topology)
{
  MetaData& meta = MetaData::get(part);

  ThrowRequireMsg(meta.is_initialized(),"set_cell_topology: initialize() must be called before this function");

  Part &root_part = meta.get_cell_topology_root_part(cell_topology);
  meta.declare_part_subset(root_part, part);
}

const std::vector<std::string>&
entity_rank_names()
{
  // TODO - Not thread safe; Use c++11 to initialize vector once c++11 is available
  static std::vector< std::string > names;
  if (names.empty()) {
    names.reserve( 4 );
    names.push_back(std::string("NODE"));
    names.push_back(std::string("EDGE"));
    names.push_back(std::string("FACE"));
    names.push_back(std::string("ELEMENT"));
    // TODO - Add constraint?
  }
  return names;
}


stk::topology
get_topology(const MetaData& meta_data, EntityRank entity_rank, const std::pair<const unsigned*, const unsigned*>& supersets)
{
  const PartVector & all_parts = meta_data.get_parts();

  stk::topology topology;

  if (supersets.first != supersets.second) {
    const Part *first_found_part = 0;

    for ( const unsigned * it = supersets.first ; it != supersets.second ; ++it ) {

      const Part & part = * all_parts[*it] ;

      if ( part.primary_entity_rank() == entity_rank ) {

        stk::topology top = meta_data.get_topology( part );

        if ( topology == stk::topology::INVALID_TOPOLOGY ) {
          topology = top ;

          if (!first_found_part)
            first_found_part = &part;
        }
        else {
          ThrowErrorMsgIf( top != stk::topology::INVALID_TOPOLOGY && top != topology,
            "topology is ambiguously defined for the bucket. It is defined as " << topology.name() <<
             " and as " << top.name() );
        }
      }
    }
  }

  return topology ;
}


stk::topology get_topology( CellTopology shards_topology, int spatial_dimension)
{
  stk::topology t;

  if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Node >()) )
    t = stk::topology::NODE;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Line<2> >()) )
    t = spatial_dimension < 2 ? stk::topology::LINE_2_1D : stk::topology::LINE_2;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Line<3> >()) )
    t = spatial_dimension < 2 ? stk::topology::LINE_3_1D : stk::topology::LINE_3;

  else if ( shards_topology == shards::getCellTopologyData< shards::Triangle<3> >() )
    t = spatial_dimension == 3 ? stk::topology::TRI_3 : stk::topology::TRI_3_2D;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()) )
    t = spatial_dimension == 3 ? stk::topology::TRI_4 : stk::topology::TRI_4_2D;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()) )
    t = spatial_dimension == 3 ? stk::topology::TRI_6 : stk::topology::TRI_6_2D;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_4 : stk::topology::QUAD_4_2D;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_8 : stk::topology::QUAD_8_2D;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_9 : stk::topology::QUAD_9_2D;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Particle >()) )
    t = stk::topology::PARTICLE;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Beam<2> >()) )
    t = stk::topology::BEAM_2;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Beam<3> >()) )
    t = stk::topology::BEAM_3;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()) )
    t = stk::topology::SHELL_LINE_2;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()) )
    t = stk::topology::SHELL_LINE_3;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()) )
    t = stk::topology::SHELL_TRI_3;
  //NOTE: shards does not define a shell triangle 4
  //else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellTriangle<4> >()) )
  //  t = stk::topology::SHELL_TRI_4;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()) )
    t = stk::topology::SHELL_TRI_6;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()) )
    t = stk::topology::SHELL_QUAD_4;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()) )
    t = stk::topology::SHELL_QUAD_8;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()) )
    t = stk::topology::SHELL_QUAD_9;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()) )
    t = stk::topology::TET_4;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()) )
    t = stk::topology::TET_8;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()) )
    t = stk::topology::TET_10;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Tetrahedron<11> >()) )
    t = stk::topology::TET_11;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()) )
    t = stk::topology::PYRAMID_5;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()) )
    t = stk::topology::PYRAMID_13;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()) )
    t = stk::topology::PYRAMID_14;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()) )
    t = stk::topology::WEDGE_6;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()) )
    t = stk::topology::WEDGE_15;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()) )
    t = stk::topology::WEDGE_18;

  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()) )
    t = stk::topology::HEX_8;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()) )
    t = stk::topology::HEX_20;
  else if ( shards_topology == CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()) )
    t = stk::topology::HEX_27;
  else if ( shards_topology.isValid() && strncmp(shards_topology.getName(), "SUPERELEMENT", 12) == 0)
    return create_superelement_topology(shards_topology.getNodeCount());

  if (t.defined_on_spatial_dimension(spatial_dimension))
    return t;

  return stk::topology::INVALID_TOPOLOGY;
}
//begin-get_cell_topology
CellTopology get_cell_topology(stk::topology t)
{
  switch(t())
  {
  case stk::topology::NODE:         
      return CellTopology( shards::getCellTopologyData< shards::Node                  >() );
  case stk::topology::LINE_2:
      return CellTopology( shards::getCellTopologyData< shards::Line<2>               >() );
  case stk::topology::LINE_3:
      return CellTopology( shards::getCellTopologyData< shards::Line<3>               >() );
  case stk::topology::TRI_3:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<3>           >() );
  case stk::topology::TRI_4:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<4>           >() );
  case stk::topology::TRI_6:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<6>           >() );
  case stk::topology::QUAD_4:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<4>      >() );
  case stk::topology::QUAD_8:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<8>      >() );
  case stk::topology::QUAD_9:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<9>      >() );
  case stk::topology::PARTICLE:
      return CellTopology( shards::getCellTopologyData< shards::Particle              >() );
  case stk::topology::LINE_2_1D:
      return CellTopology( shards::getCellTopologyData< shards::Line<2>               >() );
  case stk::topology::LINE_3_1D:
      return CellTopology( shards::getCellTopologyData< shards::Line<3>               >() );
  case stk::topology::BEAM_2:
      return CellTopology( shards::getCellTopologyData< shards::Beam<2>               >() );
  case stk::topology::BEAM_3:
      return CellTopology( shards::getCellTopologyData< shards::Beam<3>               >() );
  case stk::topology::SHELL_LINE_2:
      return CellTopology( shards::getCellTopologyData< shards::ShellLine<2>          >() );
  case stk::topology::SHELL_LINE_3:
      return CellTopology( shards::getCellTopologyData< shards::ShellLine<3>          >() );
  case stk::topology::TRI_3_2D:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<3>           >() );
  case stk::topology::TRI_4_2D:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<4>           >() );
  case stk::topology::TRI_6_2D:
      return CellTopology( shards::getCellTopologyData< shards::Triangle<6>           >() );
  case stk::topology::QUAD_4_2D:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<4>      >() );
  case stk::topology::QUAD_8_2D:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<8>      >() );
  case stk::topology::QUAD_9_2D:
      return CellTopology( shards::getCellTopologyData< shards::Quadrilateral<9>      >() );
  case stk::topology::SHELL_TRI_3:
      return CellTopology( shards::getCellTopologyData< shards::ShellTriangle<3>      >() );
  case stk::topology::SHELL_TRI_4:break;
    //NOTE: shards does not define a topology for a 4-noded triangular shell
    //return CellTopology( shards::getCellTopologyData< shards::ShellTriangle<4>    >() );
  case stk::topology::SHELL_TRI_6:
      return CellTopology( shards::getCellTopologyData< shards::ShellTriangle<6>      >() );
  case stk::topology::SHELL_QUAD_4:
      return CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<4> >() );
  case stk::topology::SHELL_QUAD_8:
      return CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<8> >() );
  case stk::topology::SHELL_QUAD_9:
      return CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<9> >() );
  case stk::topology::TET_4:
      return CellTopology( shards::getCellTopologyData< shards::Tetrahedron<4>        >() );
  case stk::topology::TET_8:
      return CellTopology( shards::getCellTopologyData< shards::Tetrahedron<8>        >() );
  case stk::topology::TET_10:
      return CellTopology( shards::getCellTopologyData< shards::Tetrahedron<10>       >() );
  case stk::topology::TET_11:
      return CellTopology( shards::getCellTopologyData< shards::Tetrahedron<11>       >() );
  case stk::topology::PYRAMID_5:
      return CellTopology( shards::getCellTopologyData< shards::Pyramid<5>            >() );
  case stk::topology::PYRAMID_13:
      return CellTopology( shards::getCellTopologyData< shards::Pyramid<13>           >() );
  case stk::topology::PYRAMID_14:
      return CellTopology( shards::getCellTopologyData< shards::Pyramid<14>           >() );
  case stk::topology::WEDGE_6:
      return CellTopology( shards::getCellTopologyData< shards::Wedge<6>              >() );
  case stk::topology::WEDGE_15:
      return CellTopology( shards::getCellTopologyData< shards::Wedge<15>             >() );
  case stk::topology::WEDGE_18:
      return CellTopology( shards::getCellTopologyData< shards::Wedge<18>             >() );
  case stk::topology::HEX_8:
      return CellTopology( shards::getCellTopologyData< shards::Hexahedron<8>         >() );
  case stk::topology::HEX_20:
      return CellTopology( shards::getCellTopologyData< shards::Hexahedron<20>        >() );
  case stk::topology::HEX_27:
      return CellTopology( shards::getCellTopologyData< shards::Hexahedron<27>        >() );
  default: break;
  }
  return CellTopology(NULL);
}
//end-get_cell_topology


FieldBase* MetaData::get_field(stk::mesh::EntityRank entity_rank, const std::string& name ) const
{
  const FieldVector& fields = m_field_repo.get_fields(static_cast<stk::topology::rank_t>(entity_rank));
  for ( FieldVector::const_iterator i =  fields.begin() ; i != fields.end(); ++i ) {
    if (equal_case((*i)->name(), name)) {
      return *i;
    }
  }
  return NULL;
}


FieldBase* get_field_by_name( const std::string& name, const MetaData & metaData )
{
  FieldBase* field = NULL;
  unsigned num_nonnull_fields = 0;
  for(stk::topology::rank_t i=stk::topology::NODE_RANK; i<=stk::topology::CONSTRAINT_RANK; ++i) {
    FieldBase* thisfield = metaData.get_field(i, name);
    if (thisfield != NULL) {
      if (field == NULL) {
        field = thisfield;
      }
      ++num_nonnull_fields;
    }
  }

  if (num_nonnull_fields > 1) {
    std::cerr << "get_field_by_name WARNING, found "<<num_nonnull_fields<<" fields with name="<<name
      <<". Returning the first one."<<std::endl;
  }

  return field;
}

void MetaData::dump_all_meta_info(std::ostream& out) const
{
  out << "MetaData info...\n";

  out << "spatial dimension = " << m_spatial_dimension << "\n";

  out << "  Entity rank names:\n";
  for (size_t i = 0, e = m_entity_rank_names.size(); i != e; ++i) {
    out << "    " << i << ": " << m_entity_rank_names[i] << "\n";
  }
  out << "  Special Parts:\n";
  out << "    Universal part ord = " << m_universal_part->mesh_meta_data_ordinal() << "\n";
  out << "    Owns part ord = " << m_owns_part->mesh_meta_data_ordinal() << "\n";
  out << "    Shared part ord = " << m_shares_part->mesh_meta_data_ordinal() << "\n";

  out << "  All parts:\n";
  const PartVector& all_parts = m_part_repo.get_all_parts();
  for(const Part* part : all_parts) {
    print(out, "    ", *part);
  }

  out << "  All fields:\n";
  const FieldVector& all_fields = m_field_repo.get_fields();
  for(const FieldBase* field : all_fields) {
     print(out, "    ", *field);
  }
}

} // namespace mesh
} // namespace stk

