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

#include <stk_mesh/base/MetaData.hpp>
#include <string.h>                     // for strcmp, strncmp
#include <Shards_CellTopologyManagedData.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <set>                          // for set
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceMin, etc
#include <stk_util/util/string_case_compare.hpp>  // for equal_case
#include "Shards_BasicTopologies.hpp"   // for getCellTopologyData, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_mesh/baseImpl/PartRepository.hpp"  // for PartRepository
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology_utils.hpp"    // for topology::num_nodes, etc
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
  for (const Part* subset : subsets) {
    if (is_topology_root_part( *subset )) {
      return true;
    }
  }
  return false;
}

void find_topologies_in_part_and_subsets_of_same_rank(const Part & part, EntityRank rank, std::set<stk::topology> & topologies_found)
{
  MetaData & meta = MetaData::get(part);
  stk::topology top = meta.get_topology(part);
  if ((top!=stk::topology::INVALID_TOPOLOGY && (part.primary_entity_rank() == rank))) {
    topologies_found.insert(top);
  }
  const PartVector & subsets = part.subsets();
  for (const Part* subset : subsets) {
    top = meta.get_topology(*subset);
    if (top!=stk::topology::INVALID_TOPOLOGY && ( (*subset).primary_entity_rank() == rank) ) {
      topologies_found.insert(top);
    }
  }
}

} // namespace

void MetaData::assign_topology(Part& part, stk::topology stkTopo)
{
  const size_t part_ordinal = part.mesh_meta_data_ordinal();

  if (part_ordinal >= m_partTopologyVector.size()) {
    m_partTopologyVector.resize(part_ordinal + 1);
  }

  m_partTopologyVector[part_ordinal] = stkTopo;

  part.m_partImpl.set_topology(stkTopo);

  ThrowRequireMsg(stkTopo != stk::topology::INVALID_TOPOLOGY, "bad topology in MetaData::assign_topology");
}

void MetaData::set_mesh_on_fields(BulkData* bulk)
{
  const FieldVector& fields = get_fields();
  for(size_t i=0; i<fields.size(); ++i) {
    fields[i]->set_mesh(bulk);
  }
}

const MetaData & MetaData::get( const BulkData & bulk_data) {
  return bulk_data.mesh_meta_data();
}

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

//----------------------------------------------------------------------

MetaData::MetaData(size_t spatial_dimension, const std::vector<std::string>& entity_rank_names)
  : m_bulk_data(NULL),
    m_commit( false ),
    m_are_late_fields_enabled( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_aura_part(NULL),
    m_field_repo(),
    m_coord_field(NULL),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/),
    m_surfaceToBlock()
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
    m_are_late_fields_enabled( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_aura_part(NULL),
    m_field_repo(),
    m_coord_field(NULL),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/),
    m_surfaceToBlock()
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_internal_part("OWNS");
  m_shares_part = & declare_internal_part("SHARES");
  m_aura_part = & declare_internal_part("AURA");
}

//----------------------------------------------------------------------

void MetaData::initialize(size_t spatial_dimension,
                          const std::vector<std::string> &rank_names,
                          const std::string &coordinate_field_name)
{
  ThrowErrorMsgIf( !m_entity_rank_names.empty(), "already initialized");
  ThrowErrorMsgIf( spatial_dimension == 0, "Min spatial dimension is 1");
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

  if (!coordinate_field_name.empty()) {
    m_coord_field_name = coordinate_field_name;
  }

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

void
MetaData::set_coordinate_field_name(const std::string & coordFieldName)
{
  m_coord_field_name = coordFieldName;
}

std::string
MetaData::coordinate_field_name() const
{
  if (!m_coord_field_name.empty()) {
    return m_coord_field_name;
  }
  else {
    return "coordinates";  // Default if nothing set by the user
  }
}

void
MetaData::set_coordinate_field(FieldBase* coord_field)
{
  m_coord_field = coord_field;
  m_coord_field_name = m_coord_field->name();
}

stk::mesh::FieldBase* try_to_find_coord_field(const stk::mesh::MetaData& meta)
{
  //attempt to initialize the coordinate-field pointer, trying a couple
  //of commonly-used names. It is expected that the client code will initialize
  //the coordinates field using set_coordinate_field, but this is an
  //attempt to be helpful for existing client codes which aren't yet calling that.

  stk::mesh::FieldBase* coord_field = meta.get_field(stk::topology::NODE_RANK, "mesh_model_coordinates");

  if (coord_field == nullptr) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "mesh_model_coordinates_0");
  }
  if (coord_field == nullptr) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "model_coordinates");
  }
  if (coord_field == nullptr) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "model_coordinates_0");
  }
  if (coord_field == nullptr) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, "coordinates");
  }
  if (coord_field == nullptr) {
    coord_field = meta.get_field(stk::topology::NODE_RANK, meta.coordinate_field_name());
  }

  return coord_field;
}

FieldBase const* MetaData::coordinate_field() const
{
  if (m_coord_field == nullptr) {
    if (!m_coord_field_name.empty()) {
      m_coord_field = get_field(stk::topology::NODE_RANK, m_coord_field_name);
    }
    else {
      m_coord_field = try_to_find_coord_field(*this);
      if (m_coord_field != nullptr) {
        m_coord_field_name = m_coord_field->name();
      }
    }
  }

  ThrowErrorMsgIf( m_coord_field == nullptr,
                   "MetaData::coordinate_field: Coordinate field has not been defined" );

  return m_coord_field;
}

//----------------------------------------------------------------------

Part * MetaData::get_part( const std::string & p_name ,
                           const char * required_by ) const
{
  Part *part = m_part_repo.get_part_by_name(p_name);
  ThrowErrorMsgIf( required_by && NULL == part,
                   "Failed to find part with name " << p_name <<
                   " for method " << required_by );
  return part;
}

Part & MetaData::declare_part( const std::string & p_name )
{
  const EntityRank rank = InvalidEntityRank;

  return *m_part_repo.declare_part( p_name, rank );
}

const char** MetaData::reserved_state_suffix() const
{
  static const char* s_reserved_state_suffix[6] = {
    "_STKFS_OLD",
    "_STKFS_N",
    "_STKFS_NM1",
    "_STKFS_NM2",
    "_STKFS_NM3",
    "_STKFS_NM4"
  };

  return s_reserved_state_suffix;
}

Part & MetaData::declare_internal_part( const std::string & p_name )
{
  std::string internal_name = impl::convert_to_internal_name(p_name);
  return declare_part(internal_name);
}

Part & MetaData::declare_part( const std::string & p_name , EntityRank rank, bool arg_force_no_induce )
{
  require_valid_entity_rank(rank);

  return *m_part_repo.declare_part( p_name , rank, arg_force_no_induce );
}

Part & MetaData::declare_internal_part( const std::string & p_name , EntityRank rank )
{
  std::string internal_name = impl::convert_to_internal_name(p_name);
  return declare_part(internal_name, rank);
}

void MetaData::declare_part_subset( Part & superset , Part & subset, bool verifyFieldRestrictions )
{
  if (!is_initialized()) {
    // can't do any topology stuff yet
    return internal_declare_part_subset(superset, subset, verifyFieldRestrictions);
  }

  stk::topology superset_stkTopo = get_topology(superset);

  const bool no_superset_topology = (superset_stkTopo == stk::topology::INVALID_TOPOLOGY);
  if ( no_superset_topology ) {
    internal_declare_part_subset(superset,subset, verifyFieldRestrictions);
    return;
  }

  // Check for topology root parts in subset or subset's subsets
  const bool subset_has_root_part = root_part_in_subset(subset);
  ThrowErrorMsgIf( subset_has_root_part, "MetaData::declare_part_subset:  Error, root cell topology part found in subset or below." );

  std::set<stk::topology> topologies;
  find_topologies_in_part_and_subsets_of_same_rank(subset,superset.primary_entity_rank(),topologies);

  ThrowErrorMsgIf( topologies.size() > 1,
      "MetaData::declare_part_subset:  Error, multiple topologies of rank "
      << superset.primary_entity_rank() << " defined below subset"
      );
  const bool non_matching_topology = ((topologies.size() == 1) && (*topologies.begin() != superset_stkTopo));
  ThrowErrorMsgIf( non_matching_topology,
      "MetaData::declare_part_subset:  Error, superset topology = "
      << superset_stkTopo.name() << " does not match the topology = "
      << topologies.begin()->name()
      << " coming from the subset part");

  // Everything is Okay!
  internal_declare_part_subset(superset,subset, verifyFieldRestrictions);
  // Update PartTopologyVector for "subset" and same-rank subsets, ad nauseum
  if (subset.primary_entity_rank() == superset.primary_entity_rank()) {
    assign_topology(subset, superset_stkTopo);
    const PartVector & subset_parts = subset.subsets();
    for (Part* part : subset_parts) {
      if (part->primary_entity_rank() == superset.primary_entity_rank()) {
        assign_topology(*part, superset_stkTopo);
      }
    }
  }
}

void MetaData::internal_declare_part_subset( Part & superset , Part & subset, bool verifyFieldRestrictions )
{
  require_not_committed();
  require_same_mesh_meta_data( MetaData::get(superset) );
  require_same_mesh_meta_data( MetaData::get(subset) );

  m_part_repo.declare_subset( superset, subset );

  if (verifyFieldRestrictions)
  {
    // The new superset / subset relationship can cause a
    // field restriction to become incompatible or redundant.
    m_field_repo.verify_and_clean_restrictions(superset, subset);
  }
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

  if (is_commit()) {
    m_bulk_data->reallocate_field_data(arg_field);
  }
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

  if (is_commit()) {
    m_bulk_data->reallocate_field_data(arg_field);
  }
}

//----------------------------------------------------------------------

void MetaData::commit()
{
  require_not_committed();

  m_commit = true ; // Cannot add or change parts or fields now

  set_mesh_on_fields(m_bulk_data);

#ifdef STK_VERBOSE_OUTPUT
  dump_all_meta_info(std::cout);
#endif
}

MetaData::~MetaData()
{
  // Destroy the properties, used 'new' to allocate so now use 'delete'

  try {
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
  // Load up appropriate standard topologies.
  register_topology(stk::topology::NODE);

  if (m_spatial_dimension == 1) {

    register_topology(stk::topology::PARTICLE);

    register_topology(stk::topology::LINE_2_1D);
    register_topology(stk::topology::LINE_3_1D);

    register_topology(stk::topology::SPRING_2);
    register_topology(stk::topology::SPRING_3);
  }

  else if (m_spatial_dimension == 2) {

    register_topology(stk::topology::LINE_2);
    register_topology(stk::topology::LINE_3);

    register_topology(stk::topology::PARTICLE);

    register_topology(stk::topology::TRI_3_2D);
    register_topology(stk::topology::TRI_4_2D);
    register_topology(stk::topology::TRI_6_2D);

    register_topology(stk::topology::QUAD_4_2D);
    register_topology(stk::topology::QUAD_8_2D);
    register_topology(stk::topology::QUAD_9_2D);

    register_topology(stk::topology::BEAM_2);
    register_topology(stk::topology::BEAM_3);

    register_topology(stk::topology::SHELL_LINE_2);
    register_topology(stk::topology::SHELL_LINE_3);

    register_topology(stk::topology::SPRING_2);
    register_topology(stk::topology::SPRING_3);
  }

  else if (m_spatial_dimension == 3) {

    register_topology(stk::topology::LINE_2);
    register_topology(stk::topology::LINE_3);

    register_topology(stk::topology::TRI_3);
    register_topology(stk::topology::TRI_4);
    register_topology(stk::topology::TRI_6);

    register_topology(stk::topology::QUAD_4);
    register_topology(stk::topology::QUAD_6);
    register_topology(stk::topology::QUAD_8);
    register_topology(stk::topology::QUAD_9);

    register_topology(stk::topology::PARTICLE);

    register_topology(stk::topology::BEAM_2);
    register_topology(stk::topology::BEAM_3);

    register_topology(stk::topology::SPRING_2);
    register_topology(stk::topology::SPRING_3);

    register_topology(stk::topology::TET_4);
    register_topology(stk::topology::TET_8);
    register_topology(stk::topology::TET_10);
    register_topology(stk::topology::TET_11);

    register_topology(stk::topology::PYRAMID_5);
    register_topology(stk::topology::PYRAMID_13);
    register_topology(stk::topology::PYRAMID_14);

    register_topology(stk::topology::WEDGE_6);
    register_topology(stk::topology::WEDGE_12);
    register_topology(stk::topology::WEDGE_15);
    register_topology(stk::topology::WEDGE_18);

    register_topology(stk::topology::HEX_8);
    register_topology(stk::topology::HEX_20);
    register_topology(stk::topology::HEX_27);

    register_topology(stk::topology::SHELL_TRI_3);
    register_topology(stk::topology::SHELL_TRI_6);

    register_topology(stk::topology::SHELL_QUAD_4);
    register_topology(stk::topology::SHELL_QUAD_8);
    register_topology(stk::topology::SHELL_QUAD_9);
  }
}

Part& MetaData::register_topology(stk::topology stkTopo)
{
  ThrowRequireMsg(is_initialized(), "MetaData::register_topology: initialize() must be called before this function");

  TopologyPartMap::iterator iter = m_topologyPartMap.find(stkTopo);
  if (iter == m_topologyPartMap.end()) {
    std::string part_name = std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + stkTopo.name();
    ThrowErrorMsgIf(get_part(part_name) != 0, "Cannot register topology with same name as existing part '" << stkTopo.name() << "'" );

    Part& part = declare_internal_part(part_name, stkTopo.rank());

    m_topologyPartMap[stkTopo] = &part;

    assign_topology(part, stkTopo);

    return part;
  }

  return *iter->second;
}

Part& MetaData::get_topology_root_part(stk::topology stkTopo) const
{
    TopologyPartMap::const_iterator iter = m_topologyPartMap.find(stkTopo);
    ThrowRequireMsg(iter != m_topologyPartMap.end(), "MetaData::get_topology_root_part ERROR, failed to map topology "<<stkTopo<<" to a part.");
    return *iter->second;
}

bool MetaData::has_topology_root_part(stk::topology stkTopo) const
{
    return (m_topologyPartMap.find(stkTopo) != m_topologyPartMap.end());
}

stk::topology MetaData::get_topology(const Part & part) const
{
  ThrowRequireMsg(is_initialized(),"MetaData::get_topology(part): initialize() must be called before this function");

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  if (part_ordinal < m_partTopologyVector.size())
  {
      return m_partTopologyVector[part_ordinal];
  }

  return stk::topology::INVALID_TOPOLOGY;
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
  stk::topology top = meta.get_topology(part);
  if (top.is_valid()) {
    const Part & root_part = meta.get_topology_root_part(top);
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
  ThrowRequireMsg(meta.is_initialized(),"set_topology: initialize() must be called before this function");

  if (part.primary_entity_rank() == InvalidEntityRank) {
    //declare_part will set the rank on the part, if the part already exists.
    meta.declare_part(part.name(), topo.rank());
  }

  Part* root_part = nullptr;
  try {
      root_part = &meta.get_topology_root_part(topo);
  }
  catch(std::exception& e) {
      meta.register_topology(topo);
      root_part = &meta.get_topology_root_part(topo);
  }

  meta.declare_part_subset(*root_part, part);
}

const std::vector<std::string>&
entity_rank_names()
{
  static std::vector< std::string > names = { std::string("NODE"), std::string("EDGE"), std::string("FACE"), std::string("ELEMENT") } ;
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
          if ( top != stk::topology::INVALID_TOPOLOGY && top != topology) {
              std::ostringstream os;
              os << "topology defined as both " << topology.name() << " and as " << top.name()
                 << "; a given mesh entity must have only one topology.";
              throw std::runtime_error(os.str());
          }
        }
      }
    }
  }

  return topology ;
}


stk::topology get_topology( shards::CellTopology shards_topology, unsigned spatial_dimension)
{
  stk::topology t;

  if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Node >()) )
    t = stk::topology::NODE;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Line<2> >()) )
    t = spatial_dimension < 2 ? stk::topology::LINE_2_1D : stk::topology::LINE_2;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Line<3> >()) )
    t = spatial_dimension < 2 ? stk::topology::LINE_3_1D : stk::topology::LINE_3;

  else if ( shards_topology == shards::getCellTopologyData< shards::Triangle<3> >() )
    t = spatial_dimension == 3 ? stk::topology::TRI_3 : stk::topology::TRI_3_2D;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()) )
    t = spatial_dimension == 3 ? stk::topology::TRI_4 : stk::topology::TRI_4_2D;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()) )
    t = spatial_dimension == 3 ? stk::topology::TRI_6 : stk::topology::TRI_6_2D;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_4 : stk::topology::QUAD_4_2D;

  //NOTE: shards does not define a quad 6
  // else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<6> >()) )
  //   t = stk::topology::QUAD_6;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_8 : stk::topology::QUAD_8_2D;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()) )
    t = spatial_dimension == 3 ? stk::topology::QUAD_9 : stk::topology::QUAD_9_2D;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Particle >()) )
    t = stk::topology::PARTICLE;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Beam<2> >()) )
    t = stk::topology::BEAM_2;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Beam<3> >()) )
    t = stk::topology::BEAM_3;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()) )
    t = stk::topology::SHELL_LINE_2;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()) )
    t = stk::topology::SHELL_LINE_3;

  //NOTE: shards does not define a spring 2
  //else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Spring<2> >()) )
  //  t = stk::topology::SPRING_2;

  //NOTE: shards does not define a spring 3
  //else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Spring<3> >()) )
  //  t = stk::topology::SPRING_3;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()) )
    t = stk::topology::SHELL_TRI_3;

  //NOTE: shards does not define a shell triangle 4
  //else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<4> >()) )
  //  t = stk::topology::SHELL_TRI_4;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()) )
    t = stk::topology::SHELL_TRI_6;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()) )
    t = stk::topology::SHELL_QUAD_4;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()) )
    t = stk::topology::SHELL_QUAD_8;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()) )
    t = stk::topology::SHELL_QUAD_9;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()) )
    t = stk::topology::TET_4;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()) )
    t = stk::topology::TET_8;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()) )
    t = stk::topology::TET_10;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<11> >()) )
    t = stk::topology::TET_11;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()) )
    t = stk::topology::PYRAMID_5;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()) )
    t = stk::topology::PYRAMID_13;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()) )
    t = stk::topology::PYRAMID_14;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()) )
    t = stk::topology::WEDGE_6;

  //NOTE: shards does not define a wedge 12
  // else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Wedge<12> >()) )
  //   t = stk::topology::WEDGE_12;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()) )
    t = stk::topology::WEDGE_15;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()) )
    t = stk::topology::WEDGE_18;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()) )
    t = stk::topology::HEX_8;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()) )
    t = stk::topology::HEX_20;
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()) )
    t = stk::topology::HEX_27;
  else if ( shards_topology.isValid() && strncmp(shards_topology.getName(), "SUPERELEMENT", 12) == 0)
    return create_superelement_topology(shards_topology.getNodeCount());
  else if ( shards_topology.isValid() && strncmp(shards_topology.getName(), "SUPERFACE", 9) == 0)
    return create_superface_topology(shards_topology.getNodeCount());
  else if ( shards_topology.isValid() && strncmp(shards_topology.getName(), "SUPEREDGE", 9) == 0)
    return create_superedge_topology(shards_topology.getNodeCount());

  if (t.defined_on_spatial_dimension(spatial_dimension))
    return t;

  return stk::topology::INVALID_TOPOLOGY;
}
//begin-get_cell_topology
shards::CellTopology get_cell_topology(stk::topology t)
{
  switch(t())
  {
  case stk::topology::NODE:         
      return shards::CellTopology( shards::getCellTopologyData< shards::Node                  >() );
  case stk::topology::LINE_2:
      return shards::CellTopology( shards::getCellTopologyData< shards::Line<2>               >() );
  case stk::topology::LINE_3:
      return shards::CellTopology( shards::getCellTopologyData< shards::Line<3>               >() );
  case stk::topology::TRI_3:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<3>           >() );
  case stk::topology::TRI_4:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<4>           >() );
  case stk::topology::TRI_6:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<6>           >() );
  case stk::topology::QUAD_4:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<4>      >() );
  case stk::topology::QUAD_6:
      //NOTE: shards does not define a topology for a 6-noded quadrilateral element
      // return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<6>      >() );
  case stk::topology::QUAD_8:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<8>      >() );
  case stk::topology::QUAD_9:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<9>      >() );
  case stk::topology::PARTICLE:
      return shards::CellTopology( shards::getCellTopologyData< shards::Particle              >() );
  case stk::topology::LINE_2_1D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Line<2>               >() );
  case stk::topology::LINE_3_1D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Line<3>               >() );
  case stk::topology::BEAM_2:
      return shards::CellTopology( shards::getCellTopologyData< shards::Beam<2>               >() );
  case stk::topology::BEAM_3:
      return shards::CellTopology( shards::getCellTopologyData< shards::Beam<3>               >() );
  case stk::topology::SHELL_LINE_2:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellLine<2>          >() );
  case stk::topology::SHELL_LINE_3:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellLine<3>          >() );
  case stk::topology::SPRING_2: break;
    //NOTE: shards does not define a topology for a 2-noded spring element
    //return shards::CellTopology( shards::getCellTopologyData< shards::Spring<2>             >() );
  case stk::topology::SPRING_3: break;
    //NOTE: shards does not define a topology for a 3-noded spring element
    //return shards::CellTopology( shards::getCellTopologyData< shards::Spring<3>             >() );
  case stk::topology::TRI_3_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<3>           >() );
  case stk::topology::TRI_4_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<4>           >() );
  case stk::topology::TRI_6_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Triangle<6>           >() );
  case stk::topology::QUAD_4_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<4>      >() );
  case stk::topology::QUAD_8_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<8>      >() );
  case stk::topology::QUAD_9_2D:
      return shards::CellTopology( shards::getCellTopologyData< shards::Quadrilateral<9>      >() );
  case stk::topology::SHELL_TRI_3:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellTriangle<3>      >() );
  case stk::topology::SHELL_TRI_4: break;
    //NOTE: shards does not define a topology for a 4-noded triangular shell
    //return shards::CellTopology( shards::getCellTopologyData< shards::ShellTriangle<4>      >() );
  case stk::topology::SHELL_TRI_6:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellTriangle<6>      >() );
  case stk::topology::SHELL_QUAD_4:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<4> >() );
  case stk::topology::SHELL_QUAD_8:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<8> >() );
  case stk::topology::SHELL_QUAD_9:
      return shards::CellTopology( shards::getCellTopologyData< shards::ShellQuadrilateral<9> >() );
  case stk::topology::TET_4:
      return shards::CellTopology( shards::getCellTopologyData< shards::Tetrahedron<4>        >() );
  case stk::topology::TET_8:
      return shards::CellTopology( shards::getCellTopologyData< shards::Tetrahedron<8>        >() );
  case stk::topology::TET_10:
      return shards::CellTopology( shards::getCellTopologyData< shards::Tetrahedron<10>       >() );
  case stk::topology::TET_11:
      return shards::CellTopology( shards::getCellTopologyData< shards::Tetrahedron<11>       >() );
  case stk::topology::PYRAMID_5:
      return shards::CellTopology( shards::getCellTopologyData< shards::Pyramid<5>            >() );
  case stk::topology::PYRAMID_13:
      return shards::CellTopology( shards::getCellTopologyData< shards::Pyramid<13>           >() );
  case stk::topology::PYRAMID_14:
      return shards::CellTopology( shards::getCellTopologyData< shards::Pyramid<14>           >() );
  case stk::topology::WEDGE_6:
      return shards::CellTopology( shards::getCellTopologyData< shards::Wedge<6>              >() );
  case stk::topology::WEDGE_12:
      //NOTE: shards does not define a topology for a 12-noded wedge
      // return shards::CellTopology( shards::getCellTopologyData< shards::Wedge<12>             >() );
  case stk::topology::WEDGE_15:
      return shards::CellTopology( shards::getCellTopologyData< shards::Wedge<15>             >() );
  case stk::topology::WEDGE_18:
      return shards::CellTopology( shards::getCellTopologyData< shards::Wedge<18>             >() );
  case stk::topology::HEX_8:
      return shards::CellTopology( shards::getCellTopologyData< shards::Hexahedron<8>         >() );
  case stk::topology::HEX_20:
      return shards::CellTopology( shards::getCellTopologyData< shards::Hexahedron<20>        >() );
  case stk::topology::HEX_27:
      return shards::CellTopology( shards::getCellTopologyData< shards::Hexahedron<27>        >() );
  default: break;
  }
  return shards::CellTopology(NULL);
}
//end-get_cell_topology


FieldBase* MetaData::get_field(stk::mesh::EntityRank entity_rank, const std::string& name ) const
{
  const FieldVector& fields = m_field_repo.get_fields(static_cast<stk::topology::rank_t>(entity_rank));
  for ( FieldBase* field : fields ) {
    if (equal_case(field->name(), name)) {
      return field;
    }
  }
  return nullptr;
}


FieldBase* get_field_by_name( const std::string& name, const MetaData & metaData )
{
  FieldBase* field = NULL;
  unsigned num_nonnull_fields = 0;
  for(stk::topology::rank_t i=stk::topology::NODE_RANK; i<=stk::topology::CONSTRAINT_RANK; ++i) {
    FieldBase* thisfield = metaData.get_field(i, name);
    if (thisfield != nullptr) {
      if (field == nullptr) {
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
  out << "MetaData info...(ptr=" << this << ")\n";

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

