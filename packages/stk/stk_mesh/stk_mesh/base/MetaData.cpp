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
#include <cstring>                      // for strcmp, strncmp
#include <Shards_CellTopologyManagedData.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <set>                          // for set
#include <algorithm>                    // for transform
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceMin, etc
#include "Shards_BasicTopologies.hpp"   // for getCellTopologyData, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Part.hpp"       // for Part, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_mesh/baseImpl/PartRepository.hpp"  // for PartRepository
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/string_case_compare.hpp>
#include <stk_mesh/base/DumpMeshInfo.hpp>

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

  part.set_topology(stkTopo);

  STK_ThrowRequireMsg(stkTopo != stk::topology::INVALID_TOPOLOGY, "bad topology in MetaData::assign_topology");
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
  STK_ThrowRequireMsg(!m_commit, "mesh MetaData has been committed.");
}

void MetaData::require_committed() const
{
  STK_ThrowRequireMsg(m_commit, "mesh MetaData has not been committed.");
}

void MetaData::require_same_mesh_meta_data( const MetaData & rhs ) const
{
  STK_ThrowRequireMsg(this == &rhs, "Different mesh_meta_data.");
}

void MetaData::require_valid_entity_rank( EntityRank rank ) const
{
  STK_ThrowRequireMsg(check_rank(rank),
      "entity_rank " << rank << " >= " << m_entity_rank_names.size() );
  STK_ThrowRequireMsg( !(rank == stk::topology::FACE_RANK && spatial_dimension() == 2),
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
    m_field_repo(*this),
    m_coord_field(NULL),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/),
    m_surfaceToBlock()
{
  const size_t numRanks = stk::topology::NUM_RANKS;
  STK_ThrowRequireMsg(entity_rank_names.size() <= numRanks, "MetaData: number of entity-ranks (" << entity_rank_names.size() << ") exceeds limit of stk::topology::NUM_RANKS (" << numRanks <<")");

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
    m_field_repo(*this),
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
  STK_ThrowErrorMsgIf( !m_entity_rank_names.empty(), "already initialized");
  STK_ThrowErrorMsgIf( spatial_dimension == 0, "Min spatial dimension is 1");
  STK_ThrowErrorMsgIf( spatial_dimension > 3, "Max spatial dimension is 3");

  if ( rank_names.empty() ) {
    m_entity_rank_names = stk::mesh::entity_rank_names();
  }
  else {
    STK_ThrowErrorMsgIf(rank_names.size() < stk::topology::ELEMENT_RANK+1,
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
  STK_ThrowErrorMsgIf( entity_rank >= entity_rank_count(),
      "entity-rank " << entity_rank <<
      " out of range. Must be in range 0.." << entity_rank_count());

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

  STK_ThrowErrorMsgIf( m_coord_field == nullptr,
                   "MetaData::coordinate_field: Coordinate field has not been defined" );

  return m_coord_field;
}

//----------------------------------------------------------------------

Part * MetaData::get_part( const std::string & p_name ,
                           const char * required_by ) const
{
  Part *part = nullptr;

  const auto iter = m_partAlias.find(p_name);

  if(iter != m_partAlias.end()) {
    part = & get_part((*iter).second);
  } else {
    part = m_part_repo.get_part_by_name(p_name);
  }

  STK_ThrowErrorMsgIf( required_by && nullptr == part,
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
  STK_ThrowRequireMsg(is_initialized(), "MetaData: Can't declare ranked part until spatial dimension has been set.");
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
  STK_ThrowErrorMsgIf( subset_has_root_part, "MetaData::declare_part_subset:  Error, root cell topology part found in subset or below." );

  std::set<stk::topology> topologies;
  find_topologies_in_part_and_subsets_of_same_rank(subset,superset.primary_entity_rank(),topologies);

  STK_ThrowErrorMsgIf( topologies.size() > 1,
      "MetaData::declare_part_subset:  Error, multiple topologies of rank "
      << superset.primary_entity_rank() << " defined below subset"
      );
  const bool non_matching_topology = ((topologies.size() == 1) && (*topologies.begin() != superset_stkTopo));
  STK_ThrowErrorMsgIf( non_matching_topology,
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
//  require_not_committed();
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
  impl::dump_all_meta_info(*this, std::cout);
#endif
}

MetaData::~MetaData()
{
  m_bulk_data = nullptr;
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

    register_topology(stk::topology::SHELL_SIDE_BEAM_2);
    register_topology(stk::topology::SHELL_SIDE_BEAM_3);

    register_topology(stk::topology::SHELL_TRI_3);
    register_topology(stk::topology::SHELL_TRI_6);

    register_topology(stk::topology::SHELL_TRI_3_ALL_FACE_SIDES);
    register_topology(stk::topology::SHELL_TRI_6_ALL_FACE_SIDES);

    register_topology(stk::topology::SHELL_QUAD_4);
    register_topology(stk::topology::SHELL_QUAD_8);
    register_topology(stk::topology::SHELL_QUAD_9);

    register_topology(stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES);
    register_topology(stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES);
    register_topology(stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES);
  }
}

Part& MetaData::register_topology(stk::topology stkTopo)
{
  STK_ThrowRequireMsg(is_initialized(), "MetaData::register_topology: initialize() must be called before this function");

  TopologyPartMap::iterator iter = m_topologyPartMap.find(stkTopo);
  if (iter == m_topologyPartMap.end()) {
    std::string part_name = std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + stkTopo.name();
    STK_ThrowErrorMsgIf(get_part(part_name) != 0, "Cannot register topology with same name as existing part '" << stkTopo.name() << "'" );

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
    STK_ThrowRequireMsg(iter != m_topologyPartMap.end(), "MetaData::get_topology_root_part ERROR, failed to map topology "<<stkTopo<<" to a part.");
    return *iter->second;
}

bool MetaData::has_topology_root_part(stk::topology stkTopo) const
{
    return (m_topologyPartMap.find(stkTopo) != m_topologyPartMap.end());
}

stk::topology MetaData::get_topology(const Part & part) const
{
  STK_ThrowRequireMsg(is_initialized(),"MetaData::get_topology(part): initialize() must be called before this function");

  PartOrdinal part_ordinal = part.mesh_meta_data_ordinal();
  if (part_ordinal < m_partTopologyVector.size())
  {
      return m_partTopologyVector[part_ordinal];
  }

  return stk::topology::INVALID_TOPOLOGY;
}

void MetaData::add_part_alias(Part& part, const std::string& alias)
{
  const auto aliasIter = m_partAlias.find(alias);
  const unsigned partOrdinal = part.mesh_meta_data_ordinal();

  if(aliasIter == m_partAlias.end()) {
    m_partAlias[alias] = partOrdinal;
  } else {
    STK_ThrowRequireMsg((*aliasIter).second == partOrdinal, "Part alias '" << alias << "' must be assigned to unique part");
  }

  std::map<unsigned, std::vector<std::string> >::iterator reverseAliasIter = m_partReverseAlias.find(partOrdinal);
  if(reverseAliasIter == m_partReverseAlias.end()) {
    std::vector<std::string> entry{alias};
    m_partReverseAlias[partOrdinal] = entry;
  }
  else {
    stk::util::insert_keep_sorted_and_unique(alias, (*reverseAliasIter).second);
  }
}

bool MetaData::delete_part_alias(Part& part, const std::string& alias)
{
  bool deleted = false;
  unsigned partOrdinal = stk::mesh::InvalidOrdinal;
  auto iter = m_partAlias.find(alias);

  if(iter != m_partAlias.end()) {
    partOrdinal = iter->second;
    m_partAlias.erase(iter);

    std::map<unsigned, std::vector<std::string> >::iterator reverseAliasIter = m_partReverseAlias.find(partOrdinal);
    STK_ThrowRequireMsg(reverseAliasIter != m_partReverseAlias.end(), "Could not find reverse alias map entry for part: " << part.name());

    std::vector<std::string>& aliases = reverseAliasIter->second;
    auto result = std::lower_bound(aliases.begin(), aliases.end(), alias );
    STK_ThrowRequireMsg((result != aliases.end()) && (*result == alias),
                    "Could not find alias: '" << alias << "' for part: " << part.name() << " in reverse alias map");

    aliases.erase(result);
    deleted = true;
  }

  return deleted;
}

bool MetaData::delete_part_alias_case_insensitive(Part& part, const std::string& alias)
{
  bool deleted = false;
  unsigned partOrdinal = stk::mesh::InvalidOrdinal;
  for(auto iter = m_partAlias.begin(); iter != m_partAlias.end(); ) {
    if(stk::equal_case(iter->first, alias)) {
      STK_ThrowRequireMsg((partOrdinal == stk::mesh::InvalidOrdinal) || (partOrdinal == iter->second),
                      "Part alias '" << alias << "' not  uniquely assigned");
      partOrdinal = iter->second;
      iter = m_partAlias.erase(iter);
      deleted = true;
    } else {
      iter++;
    }
  }

  if(partOrdinal != stk::mesh::InvalidOrdinal) {
    std::map<unsigned, std::vector<std::string> >::iterator reverseAliasIter = m_partReverseAlias.find(partOrdinal);
    STK_ThrowRequireMsg(reverseAliasIter != m_partReverseAlias.end(), "Could not find reverse alias map entry for part: " << part.name());

    std::vector<std::string>& aliases = reverseAliasIter->second;
    for(auto iter = aliases.begin(); iter != aliases.end(); ) {
      if(stk::equal_case(*iter, alias)) {
        iter = aliases.erase(iter);
        deleted = true;
      } else {
        iter++;
      }
    }
  }

  return deleted;
}

std::vector<std::string> MetaData::get_part_aliases(const Part& part) const
{
  auto iter = m_partReverseAlias.find(part.mesh_meta_data_ordinal());

  if(iter != m_partReverseAlias.end())
    return (*iter).second;

  return std::vector<std::string>();
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Verify parallel consistency of fields and parts

namespace {

class VerifyConsistency
{
public:
  VerifyConsistency(const MetaData & meta, ParallelMachine pm)
    : m_meta(meta),
      m_parallelMachine(pm),
      m_pRank(stk::parallel_machine_rank(pm))
  {
  }

  virtual ~VerifyConsistency() = default;

  void verify()
  {
    const int rootRank = 0;

    CommBroadcast comm(m_parallelMachine, rootRank);
    stk::pack_and_communicate(comm, [&](){
      if (m_pRank == rootRank) {
        pack(comm.send_buffer());
      }
    });

    m_errStream.str(""); //clear
    int ok = unpack_verify(comm.recv_buffer());
    stk::all_reduce(m_parallelMachine, ReduceMin<1>(&ok));

    STK_ThrowRequireMsg(ok, "[p" << m_pRank << "] MetaData parallel consistency failure: "<<m_errStream.str());
  }

  virtual void pack(CommBuffer & b) = 0;
  virtual bool unpack_verify(CommBuffer & b) = 0;

protected:
  const MetaData & m_meta;
  ParallelMachine m_parallelMachine;
  int m_pRank;
  std::ostringstream m_errStream;
};

class VerifyPartConsistency : public VerifyConsistency
{
public:
  VerifyPartConsistency(const MetaData & meta, ParallelMachine pm)
    : VerifyConsistency(meta, pm),
      m_rootPartOrdinal(0),
      m_rootPartNameLen(0),
      m_rootPartName{0},
      m_rootPartRank(stk::topology::INVALID_RANK),
      m_rootPartTopology(stk::topology::INVALID_TOPOLOGY),
      m_rootPartSubsetSize(0)
  {
  }

  virtual ~VerifyPartConsistency() = default;

  virtual void pack(CommBuffer & b)
  {
    for (const stk::mesh::Part * part : m_meta.get_parts()) {
      if (!stk::mesh::impl::is_internal_part(*part)) {
        const PartVector & subsetParts = part->subsets();
        const char * const partNamePtr = part->name().c_str();
        const unsigned     partNameLen = static_cast<unsigned>(part->name().size()) + 1;

        b.pack<stk::mesh::Ordinal>(part->mesh_meta_data_ordinal());
        b.pack<unsigned>(partNameLen);
        b.pack<char>(partNamePtr, partNameLen);
        b.pack<stk::mesh::EntityRank>(part->primary_entity_rank());
        b.pack<stk::topology::topology_t>(part->topology());
        b.pack<unsigned>(subsetParts.size());
        for (const stk::mesh::Part * subsetPart : subsetParts) {
          b.pack<unsigned>(subsetPart->mesh_meta_data_ordinal());
        }
      }
    }
  }

  virtual bool unpack_verify(CommBuffer & b)
  {
    const stk::mesh::PartVector & localParts = m_meta.get_parts();
    bool ok = true;

    std::set<const stk::mesh::Part *> unprocessedParts;
    std::transform(localParts.begin(), localParts.end(),
                   std::inserter(unprocessedParts, unprocessedParts.end()), [](const stk::mesh::Part* p) {return p;});

    while (b.remaining()) {
      b.unpack<stk::mesh::Ordinal>(m_rootPartOrdinal);
      b.unpack<unsigned>(m_rootPartNameLen);
      b.unpack<char>(m_rootPartName, m_rootPartNameLen);
      b.unpack<stk::mesh::EntityRank>(m_rootPartRank);
      b.unpack<stk::topology::topology_t>(m_rootPartTopology);
      b.unpack<unsigned>(m_rootPartSubsetSize);
      m_rootPartSubsetOrdinals.resize(m_rootPartSubsetSize);
      for (unsigned i = 0; i < m_rootPartSubsetSize; ++i) {
        b.unpack<unsigned>(m_rootPartSubsetOrdinals[i]);
      }

      stk::mesh::Part * localPart = get_local_part(localParts);
      unprocessedParts.erase(localPart);
      ok = check_local_part(localPart) && ok;

      ok = ok && check_part_name(localPart);
      ok = ok && check_part_rank(localPart);
      ok = ok && check_part_topology(localPart);
      ok = ok && check_part_subsets(localPart);
    }

    ok = check_extra_parts(unprocessedParts) && ok;

    return ok;
  }

  stk::mesh::Part * get_local_part(const stk::mesh::PartVector & parts)
  {
    stk::mesh::Part * localPart = nullptr;
    if (m_rootPartOrdinal < parts.size()) {
      localPart = parts[m_rootPartOrdinal];
    }
    return localPart;
  }

  bool check_local_part(const stk::mesh::Part * localPart)
  {
    bool localPartValid = true;
    if (localPart == nullptr) {
      m_errStream << "[p" << m_pRank << "] Received extra Part (" << m_rootPartName << ") from root processor" << std::endl;
      localPartValid = false;
    }
    return localPartValid;
  }

  bool check_part_name(const stk::mesh::Part * localPart)
  {
    const char * const localPartNamePtr = localPart->name().c_str();
    const unsigned localPartNameLen = static_cast<unsigned>(localPart->name().size()) + 1;
    bool nameLengthMatches = (localPartNameLen == m_rootPartNameLen);
    bool nameMatches = nameLengthMatches && (strncmp(localPartNamePtr, m_rootPartName, localPartNameLen-1) == 0);

    if (!nameMatches) {
      m_errStream << "[p" << m_pRank << "] Part name (" << localPart->name()
                << ") does not match Part name (" << m_rootPartName << ") on root processor" << std::endl;
    }
    return nameMatches;
  }

  bool check_part_rank(const stk::mesh::Part * localPart)
  {
    stk::mesh::EntityRank localPartRank = localPart->primary_entity_rank();
    bool partRankMatches = (localPartRank == m_rootPartRank);

    if (!partRankMatches) {
      m_errStream << "[p" << m_pRank << "] Part " << localPart->name() << " rank (" << localPartRank
                << ") does not match Part " << m_rootPartName << " rank (" << m_rootPartRank << ") on root processor" << std::endl;
    }
    return partRankMatches;
  }

  bool check_part_topology(const stk::mesh::Part * localPart)
  {
    stk::topology localPartTopology = localPart->topology();
    bool partTopologyMatches = (localPartTopology == m_rootPartTopology);

    if (!partTopologyMatches) {
      m_errStream << "[p" << m_pRank << "] Part " << localPart->name() << " topology (" << localPartTopology
                << ") does not match Part " << m_rootPartName << " topology (" << stk::topology(m_rootPartTopology)
                << ") on root processor" << std::endl;
    }
    return partTopologyMatches;
  }

  bool check_part_subsets(const stk::mesh::Part * localPart)
  {
    const PartVector & localSubsetParts = localPart->subsets();
    bool numberOfSubsetsMatches = (localSubsetParts.size() == m_rootPartSubsetSize);

    bool subsetMatches = numberOfSubsetsMatches;
    if (numberOfSubsetsMatches) {
      for (unsigned i = 0; i < m_rootPartSubsetSize; ++i) {
        subsetMatches = subsetMatches && (localSubsetParts[i]->mesh_meta_data_ordinal() == m_rootPartSubsetOrdinals[i]);
      }
    }

    if (!subsetMatches) {
      m_errStream << "[p" << m_pRank << "] Part " << localPart->name() << " subset ordinals (";
      for (const stk::mesh::Part * subsetPart : localSubsetParts) {
        m_errStream << subsetPart->mesh_meta_data_ordinal() << " ";
      }
      m_errStream << ") does not match Part " << m_rootPartName << " subset ordinals (";
      for (const unsigned rootPartSubsetOrdinal : m_rootPartSubsetOrdinals) {
        m_errStream << rootPartSubsetOrdinal << " ";
      }
      m_errStream << ") on root processor" << std::endl;
    }

    return subsetMatches;
  }

  bool check_extra_parts(const std::set<const stk::mesh::Part *> & unprocessedParts)
  {
    bool noExtraParts = true;
    for (const stk::mesh::Part * part : unprocessedParts) {
      if (!stk::mesh::impl::is_internal_part(*part)) {
        m_errStream << "[p" << m_pRank << "] Have extra Part (" << part->name() << ") that does not exist on root processor" << std::endl;
        noExtraParts = false;
      }
    }
    return noExtraParts;
  }

private:
  stk::mesh::Ordinal m_rootPartOrdinal;
  unsigned m_rootPartNameLen;
  char m_rootPartName[4096];
  stk::mesh::EntityRank m_rootPartRank;
  stk::topology::topology_t m_rootPartTopology;
  unsigned m_rootPartSubsetSize;
  std::vector<unsigned> m_rootPartSubsetOrdinals;
};


class VerifyFieldConsistency : public VerifyConsistency
{
public:
  VerifyFieldConsistency(const MetaData & meta, ParallelMachine pm)
    : VerifyConsistency(meta, pm),
      m_rootFieldOrdinal(0),
      m_rootFieldNameLen(0),
      m_rootFieldName{0},
      m_rootFieldRank(stk::topology::INVALID_RANK),
      m_rootFieldNumberOfStates(0)
  {
  }

  virtual ~VerifyFieldConsistency() = default;

  virtual void pack(CommBuffer & b)
  {
    for (const stk::mesh::FieldBase * field : m_meta.get_fields()) {
      const char * const fieldNamePtr = field->name().c_str();
      const unsigned     fieldNameLen = static_cast<unsigned>(field->name().size()) + 1;

      b.pack<stk::mesh::Ordinal>(field->mesh_meta_data_ordinal());
      b.pack<unsigned>(fieldNameLen);
      b.pack<char>(fieldNamePtr, fieldNameLen);
      b.pack<stk::mesh::EntityRank>(field->entity_rank());
      b.pack<unsigned>(field->number_of_states());

      //  Remaining Fields attributes that should be checked:
      //    field->data_traits()
      //    field->state()
      //    field->restrictions()
      //    field->attribute()
    }
  }

  virtual bool unpack_verify(CommBuffer & b)
  {
    const std::vector<stk::mesh::FieldBase *> & localFields = m_meta.get_fields();
    bool ok = true;

    std::set<const stk::mesh::FieldBase *> unprocessedFields;
    std::transform(localFields.begin(), localFields.end(),
                   std::inserter(unprocessedFields, unprocessedFields.end()), [](const stk::mesh::FieldBase* p) {return p;});

    while (b.remaining()) {
      b.unpack<stk::mesh::Ordinal>(m_rootFieldOrdinal);
      b.unpack<unsigned>(m_rootFieldNameLen);
      b.unpack<char>(m_rootFieldName, m_rootFieldNameLen);
      b.unpack<stk::mesh::EntityRank>(m_rootFieldRank);
      b.unpack<unsigned>(m_rootFieldNumberOfStates);

      const stk::mesh::FieldBase * localField = get_local_field(localFields);
      unprocessedFields.erase(localField);
      ok = check_local_field(localField) && ok;

      ok = ok && check_field_name(localField);
      ok = ok && check_field_rank(localField);
      ok = ok && check_field_number_of_states(localField);
    }

    ok = check_extra_fields(unprocessedFields) && ok;

    return ok;
  }

  stk::mesh::FieldBase * get_local_field(const std::vector<FieldBase *> & localFields)
  {
    stk::mesh::FieldBase * localField = nullptr;
    if (m_rootFieldOrdinal < localFields.size()) {
      localField = localFields[m_rootFieldOrdinal];
    }
    return localField;
  }

  bool check_local_field(const stk::mesh::FieldBase * localField)
  {
    bool localFieldValid = true;
    if (localField == nullptr) {
      m_errStream << "[p" << m_pRank << "] Received extra Field (" << m_rootFieldName << ") from root processor" << std::endl;
      localFieldValid = false;
    }
    return localFieldValid;
  }

  bool check_field_name(const stk::mesh::FieldBase * localField)
  {
    const char * const localFieldNamePtr = localField->name().c_str();
    const unsigned localFieldNameLen = static_cast<unsigned>(localField->name().size()) + 1;
    bool fieldNameLengthMatches = (localFieldNameLen == m_rootFieldNameLen);
    bool fieldNameMatches = fieldNameLengthMatches && (strncmp(localFieldNamePtr, m_rootFieldName, localFieldNameLen-1) == 0);

    if (!fieldNameMatches) {
      m_errStream << "[p" << m_pRank << "] Field name (" << localField->name()
                << ") does not match Field name (" << m_rootFieldName << ") on root processor" << std::endl;
    }
    return fieldNameMatches;
  }

  bool check_field_rank(const stk::mesh::FieldBase * localField)
  {
    stk::mesh::EntityRank localFieldRank = localField->entity_rank();
    bool fieldRankMatches = (localFieldRank == m_rootFieldRank);

    if (!fieldRankMatches) {
      m_errStream << "[p" << m_pRank << "] Field " << localField->name() << " rank (" << localFieldRank
                << ") does not match Field " << m_rootFieldName << " rank (" << m_rootFieldRank << ") on root processor" << std::endl;
    }
    return fieldRankMatches;
  }

  bool check_field_number_of_states(const stk::mesh::FieldBase * localField)
  {
    unsigned localNumberOfStates = localField->number_of_states();
    bool fieldNumberOfStatesMatches = (localNumberOfStates == m_rootFieldNumberOfStates);

    if (!fieldNumberOfStatesMatches) {
      m_errStream << "[p" << m_pRank << "] Field " << localField->name() << " number of states (" << localNumberOfStates
                << ") does not match Field " << m_rootFieldName << " number of states (" << m_rootFieldNumberOfStates
                << ") on root processor" << std::endl;
    }
    return fieldNumberOfStatesMatches;
  }

  bool check_extra_fields(const std::set<const stk::mesh::FieldBase *> & unprocessedFields)
  {
    bool noExtraFields = true;
    for (const stk::mesh::FieldBase * field : unprocessedFields) {
      m_errStream << "[p" << m_pRank << "] Have extra Field (" << field->name() << ") that does not exist on root processor" << std::endl;
      noExtraFields = false;
    }
    return noExtraFields;
  }

private:
  stk::mesh::Ordinal m_rootFieldOrdinal;
  unsigned m_rootFieldNameLen;
  char m_rootFieldName[4096];
  stk::mesh::EntityRank m_rootFieldRank;
  unsigned m_rootFieldNumberOfStates;
};

}

void verify_parallel_consistency(const MetaData & meta, ParallelMachine pm)
{
  VerifyPartConsistency partConsistency(meta, pm);
  partConsistency.verify();

  VerifyFieldConsistency fieldConsistency(meta, pm);
  fieldConsistency.verify();
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
  STK_ThrowRequireMsg(meta.is_initialized(),"set_topology: initialize() must be called before this function");

  if (part.primary_entity_rank() == InvalidEntityRank) {
    //declare_part will set the rank on the part, if the part already exists.
    meta.declare_part(part.name(), topo.rank());
  }

  Part* root_part = nullptr;
  try {
      root_part = &meta.get_topology_root_part(topo);
  }
  catch(std::exception&) {
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
          STK_ThrowRequireMsg(top == stk::topology::INVALID_TOPOLOGY || top == topology,
              "topology defined as both " << topology.name() << " and as " << top.name()
                  << "; a given mesh entity must have only one topology.");
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

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()) ) {
    t = stk::topology::SHELL_TRI_3;
    // t = stk::topology::SHELL_TRI_3_ALL_FACE_SIDES;
  }

  //NOTE: shards does not define a shell triangle 4
  //else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<4> >()) )
  //  t = stk::topology::SHELL_TRI_4;

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()) ) {
    t = stk::topology::SHELL_TRI_6;
    // t = stk::topology::SHELL_TRI_6_ALL_FACE_SIDES;
  }

  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()) ) {
    t = stk::topology::SHELL_QUAD_4;
  //   t = stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES;
  }
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()) ) {
    t = stk::topology::SHELL_QUAD_8;
  //   t = stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES;
  }
  else if ( shards_topology == shards::CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()) ) {
    t = stk::topology::SHELL_QUAD_9;
  //   t = stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES;
  }

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
    return shards::CellTopology(shards::getCellTopologyData<shards::Node>());
  case stk::topology::LINE_2:
    return shards::CellTopology(shards::getCellTopologyData<shards::Line<2>>());
  case stk::topology::LINE_3:
    return shards::CellTopology(shards::getCellTopologyData<shards::Line<3>>());
  case stk::topology::TRI_3:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3>>());
  case stk::topology::TRI_4:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<4>>());
  case stk::topology::TRI_6:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<6>>());
  case stk::topology::QUAD_4:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>());
  case stk::topology::QUAD_6: break;
    //NOTE: shards does not define a topology for a 6-noded quadrilateral element
    // return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<6>>());
  case stk::topology::QUAD_8:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<8>>());
  case stk::topology::QUAD_9:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<9>>());
  case stk::topology::PARTICLE:
    return shards::CellTopology(shards::getCellTopologyData<shards::Particle>());
  case stk::topology::LINE_2_1D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Line<2>>());
  case stk::topology::LINE_3_1D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Line<3>>());
  case stk::topology::BEAM_2:
    return shards::CellTopology(shards::getCellTopologyData<shards::Beam<2>>());
  case stk::topology::BEAM_3:
    return shards::CellTopology(shards::getCellTopologyData<shards::Beam<3>>());
  case stk::topology::SHELL_LINE_2:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<2>>());
  case stk::topology::SHELL_LINE_3:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<3>>());
  case stk::topology::SPRING_2: break;
    //NOTE: shards does not define a topology for a 2-noded spring element
    //return shards::CellTopology(shards::getCellTopologyData<shards::Spring<2>>());
  case stk::topology::SPRING_3: break;
    //NOTE: shards does not define a topology for a 3-noded spring element
    //return shards::CellTopology(shards::getCellTopologyData<shards::Spring<3>>());
  case stk::topology::TRI_3_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3>>());
  case stk::topology::TRI_4_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<4>>());
  case stk::topology::TRI_6_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Triangle<6>>());
  case stk::topology::QUAD_4_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>());
  case stk::topology::QUAD_8_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<8>>());
  case stk::topology::QUAD_9_2D:
    return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<9>>());
  case stk::topology::SHELL_TRI_3:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<3>>());
  case stk::topology::SHELL_TRI_4: break;
    //NOTE: shards does not define a topology for a 4-noded triangular shell
    //return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<4>>());
  case stk::topology::SHELL_TRI_6:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<6>>());
  case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<3>>());
  case stk::topology::SHELL_TRI_4_ALL_FACE_SIDES: break;
    //NOTE: shards does not define a topology for a 4-noded triangular shell
    //return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<4>>());
  case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<6>>());
  case stk::topology::SHELL_QUAD_4:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<4>>());
  case stk::topology::SHELL_QUAD_8:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<8>>());
  case stk::topology::SHELL_QUAD_9:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<9>>());
  case stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<4>>());
  case stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<8>>());
  case stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES:
    return shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<9>>());
  case stk::topology::TET_4:
    return shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4>>());
  case stk::topology::TET_8:
    return shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<8>>());
  case stk::topology::TET_10:
    return shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<10>>());
  case stk::topology::TET_11:
    return shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<11>>());
  case stk::topology::PYRAMID_5:
    return shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5>>());
  case stk::topology::PYRAMID_13:
    return shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<13>>());
  case stk::topology::PYRAMID_14:
    return shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<14>>());
  case stk::topology::WEDGE_6:
    return shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6>>());
  case stk::topology::WEDGE_12: break;
    //NOTE: shards does not define a topology for a 12-noded wedge
    // return shards::CellTopology(shards::getCellTopologyData<shards::Wedge<12>>());
  case stk::topology::WEDGE_15:
    return shards::CellTopology(shards::getCellTopologyData<shards::Wedge<15>>());
  case stk::topology::WEDGE_18:
    return shards::CellTopology(shards::getCellTopologyData<shards::Wedge<18>>());
  case stk::topology::HEX_8:
    return shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8>>());
  case stk::topology::HEX_20:
    return shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<20>>());
  case stk::topology::HEX_27:
    return shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<27>>());
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

void MetaData::use_simple_fields() { }

bool MetaData::is_using_simple_fields() const { return true; }

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

void sync_to_host_and_mark_modified(const MetaData& meta)
{
  const std::vector<FieldBase*>& fields = meta.get_fields();
  for(FieldBase* field : fields) {
    if (field->number_of_states() > 1) {
      field->sync_to_host();
      field->modify_on_host();
    }
  }
}

} // namespace mesh
} // namespace stk

