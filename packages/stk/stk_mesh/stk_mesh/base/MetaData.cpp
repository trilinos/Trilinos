/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <string.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CellTopology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_mesh/baseImpl/FieldRepository.hpp>

namespace stk {
namespace mesh {

namespace {

void assign_cell_topology(
  MetaData::PartCellTopologyVector & part_cell_topology_vector,
  size_t                   part_ordinal,
  const CellTopology       cell_topology)
{
  TraceIfWatching("stk::mesh::assign_cell_topology", LOG_PART, part_ordinal);
  DiagIfWatching(LOG_PART, part_ordinal, "assigning cell topo: " << cell_topology.getName());

  if (part_ordinal >= part_cell_topology_vector.size()) {
    part_cell_topology_vector.resize(part_ordinal + 1);
  }

  part_cell_topology_vector[part_ordinal] = cell_topology;

  ThrowRequireMsg(cell_topology.getCellTopologyData(), "bad topology in MetaData::assign_cell_topology");
}

bool root_part_in_subset(stk::mesh::Part & part)
{
  if (is_cell_topology_root_part(part)) {
    return true;
  }
  const PartVector & subsets = part.subsets();
  for (PartVector::const_iterator it=subsets.begin() ; it != subsets.end() ; ++it) {
    if (is_cell_topology_root_part( **it )) {
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

} // namespace

MetaData & MetaData::get( const BulkData & bulk_data) {
  return bulk_data.meta_data();
}

MetaData & MetaData::get( const Bucket & bucket) {
  return MetaData::get(BulkData::get(bucket));
}

MetaData & MetaData::get( const Entity & entity) {
  return MetaData::get(BulkData::get(entity));
}

MetaData & MetaData::get( const Ghosting & ghost) {
  return MetaData::get(BulkData::get(ghost));
}
//----------------------------------------------------------------------

std::ostream &
print_entity_id( std::ostream & os , const MetaData & meta_data ,
                  unsigned type , EntityId id )
{
  const std::string & name = meta_data.entity_rank_name( type );
  return os << name << "[" << id << "]" ;
}


std::ostream &
print_entity_key( std::ostream & os , const MetaData & meta_data ,
                  const EntityKey & key )
{
  const unsigned type   = entity_rank(key);
  const EntityId id = entity_id(key);
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
}

void MetaData::require_not_relation_target( const Part * const part ) const
{
  std::vector<PartRelation>::const_iterator i_end = part->relations().end();
  std::vector<PartRelation>::const_iterator i     = part->relations().begin();
  for ( ; i != i_end ; ++i ) {
    ThrowRequireMsg( part != i->m_target,
                     "Part[" << part->name() << "] is a PartRelation target");
  }
}

//----------------------------------------------------------------------

MetaData::MetaData(size_t spatial_dimension, const std::vector<std::string>& entity_rank_names)
  : m_commit( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_field_repo(),
    m_field_relations( ),
    m_properties( ),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/)
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_internal_part("OWNS");
  m_shares_part = & declare_internal_part("SHARES");

  initialize(spatial_dimension, entity_rank_names);
}

MetaData::MetaData()
  : m_commit( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_field_repo(),
    m_field_relations( ),
    m_properties( ),
    m_entity_rank_names( ),
    m_spatial_dimension( 0 /*invalid spatial dimension*/)
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_internal_part("OWNS");
  m_shares_part = & declare_internal_part("SHARES");
}

//----------------------------------------------------------------------

void MetaData::initialize(size_t spatial_dimension, const std::vector<std::string> &rank_names)
{
  ThrowErrorMsgIf( !m_entity_rank_names.empty(), "already initialized");

  if ( rank_names.empty() ) {
    m_entity_rank_names = stk::mesh::entity_rank_names(spatial_dimension);
  }
  else {
    ThrowErrorMsgIf(rank_names.size() < spatial_dimension+1,
                    "Entity rank name vector must name every rank, rank_names.size() = " <<
                    rank_names.size() << ", need " << spatial_dimension+1 << " names");
    m_entity_rank_names = rank_names;
  }

  m_spatial_dimension = spatial_dimension;

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
      entity_rank = i;
      break;
    }
  return entity_rank;
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

Part & MetaData::declare_part( const std::string & p_name )
{
  require_not_committed();

  const EntityRank rank = InvalidEntityRank;

  return *m_part_repo.declare_part( p_name, rank );
}

Part & MetaData::declare_internal_part( const std::string & p_name )
{
  std::string internal_name = convert_to_internal_name(p_name);
  return declare_part(internal_name);
}

Part & MetaData::declare_part( const std::string & p_name , EntityRank rank )
{
  require_not_committed();
  require_valid_entity_rank(rank);

  return *m_part_repo.declare_part( p_name , rank );
}

Part & MetaData::declare_internal_part( const std::string & p_name , EntityRank rank )
{
  std::string internal_name = convert_to_internal_name(p_name);
  return declare_part(internal_name, rank);
}

Part & MetaData::declare_part( const PartVector & part_intersect )
{
  require_not_committed();

  for ( PartVector::const_iterator
        i = part_intersect.begin() ; i != part_intersect.end() ; ++i ) {
    require_not_relation_target(*i);
  }

  return *m_part_repo.declare_part( part_intersect );
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
    assign_cell_topology(m_partCellTopologyVector, subset.mesh_meta_data_ordinal(), superset_top);
    const PartVector & subset_parts = subset.subsets();
    for (PartVector::const_iterator it=subset_parts.begin() ; it != subset_parts.end() ; ++it) {
      const Part & it_part = **it;
      if (it_part.primary_entity_rank() == superset.primary_entity_rank()) {
        assign_cell_topology(m_partCellTopologyVector, it_part.mesh_meta_data_ordinal(), superset_top);
      }
    }
  }
}

void MetaData::internal_declare_part_subset( Part & superset , Part & subset )
{
  static const char method[] = "stk::mesh::MetaData::declare_part_subset" ;

  require_not_committed();
  require_same_mesh_meta_data( MetaData::get(superset) );
  require_same_mesh_meta_data( MetaData::get(subset) );
  require_not_relation_target( &superset );
  require_not_relation_target( &subset );

  m_part_repo.declare_subset( superset, subset );

  // The new superset / subset relationship can cause a
  // field restriction to become incompatible or redundant.
  m_field_repo.verify_and_clean_restrictions(method, superset, subset, m_part_repo.get_all_parts());
}

void MetaData::declare_part_relation(
  Part & root_part ,
  relation_stencil_ptr stencil ,
  Part & target_part )
{
  require_not_committed();
  require_not_relation_target( &root_part );

  ThrowErrorMsgIf( !stencil, "stencil function pointer cannot be NULL" );

  ThrowErrorMsgIf( 0 != target_part.subsets().size() ||
                   0 != target_part.intersection_of().size() ||
                   1 != target_part.supersets().size(),
                   "target Part[" << target_part.name() <<
                   "] cannot be a superset or subset" );

  PartRelation tmp ;
  tmp.m_root = & root_part ;
  tmp.m_target = & target_part ;
  tmp.m_function = stencil ;

  m_part_repo.declare_part_relation( root_part, tmp, target_part );
}

//----------------------------------------------------------------------

FieldBase *
MetaData::declare_field_base(
  const std::string & arg_name ,
  const DataTraits  & arg_traits ,
  unsigned            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned            arg_num_states )
{
  require_not_committed();

  return m_field_repo.declare_field(
                arg_name,
                arg_traits,
                arg_rank,
                arg_dim_tags,
                arg_num_states,
                this
               );
}

void MetaData::declare_field_restriction(
  FieldBase      & arg_field ,
  EntityRank       arg_entity_rank ,
  const Part     & arg_part ,
  const unsigned * arg_stride ,
  const void     * arg_init_value )
{
  static const char method[] =
    "std::mesh::MetaData::declare_field_restriction" ;

  //require_not_committed(); // Moved to FieldBaseImpl::declare_field_restriction
  require_same_mesh_meta_data( MetaData::get(arg_field) );
  require_same_mesh_meta_data( MetaData::get(arg_part) );

  m_field_repo.declare_field_restriction(
      method,
      arg_field,
      arg_entity_rank,
      arg_part,
      m_part_repo.get_all_parts(),
      arg_stride,
      arg_init_value
      );
}

void MetaData::declare_field_restriction(
  FieldBase      & arg_field ,
  EntityRank       arg_entity_rank ,
  const Selector & arg_selector ,
  const unsigned * arg_stride ,
  const void     * arg_init_value )
{
  static const char method[] =
    "std::mesh::MetaData::declare_field_restriction" ;

  //require_not_committed(); // Moved to FieldBaseImpl::declare_field_restriction
  require_same_mesh_meta_data( MetaData::get(arg_field) );

  m_field_repo.declare_field_restriction(
      method,
      arg_field,
      arg_entity_rank,
      arg_selector,
      m_part_repo.get_all_parts(),
      arg_stride,
      arg_init_value
      );
}


void MetaData::internal_declare_field_relation(
  FieldBase & pointer_field ,
  relation_stencil_ptr stencil ,
  FieldBase & referenced_field )
{
  FieldRelation tmp ;
  tmp.m_root   = & pointer_field ;
  tmp.m_target = & referenced_field ;
  tmp.m_function = stencil ;

  m_field_relations.push_back( tmp );
}

//----------------------------------------------------------------------

void MetaData::commit()
{
  require_not_committed();

  m_commit = true ; // Cannot add or change parts or fields now
}

MetaData::~MetaData()
{
  // Destroy the properties, used 'new' to allocate so now use 'delete'

  try {
    std::vector<PropertyBase * >::iterator j = m_properties.begin();

    for ( ; j != m_properties.end() ; ++j ) { delete *j ; }

    m_properties.clear();
  } catch(...) {}

  // PartRepository is member data
  // FieldRepository is member data
}

void MetaData::internal_declare_known_cell_topology_parts()
{
  // Load up appropriate standard cell topologies.
  register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Node >()), NODE_RANK);

  if (m_spatial_dimension == 1) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), element_rank()); // ???
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), element_rank()); // ???

  }

  else if (m_spatial_dimension == 2) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellLine<2> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellLine<3> >()), element_rank());
  }

  else if (m_spatial_dimension == 3) {

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<2> >()), EDGE_RANK);
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Line<3> >()), EDGE_RANK);

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<6> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Triangle<4> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<8> >()), side_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Quadrilateral<9> >()), side_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Particle >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<2> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Beam<3> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<10> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<11> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Tetrahedron<8> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<5> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<13> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Pyramid<14> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<6> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<15> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Wedge<18> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::Hexahedron<27> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellTriangle<6> >()), element_rank());

    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<4> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<8> >()), element_rank());
    register_cell_topology(CellTopology(shards::getCellTopologyData< shards::ShellQuadrilateral<9> >()), element_rank());
  }
}

void MetaData::register_cell_topology(const CellTopology cell_topology, EntityRank entity_rank)
{
  ThrowRequireMsg(is_initialized(),"MetaData::register_cell_topology: initialize() must be called before this function");

  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);

  const bool       duplicate     = it != m_cellTopologyPartEntityRankMap.end();
  const EntityRank existing_rank = duplicate ? (*it).second.second : 0;

  ThrowInvalidArgMsgIf(m_spatial_dimension < entity_rank,
    "entity_rank " << entity_rank << ", " <<
    "exceeds maximum spatial_dimension = " << m_spatial_dimension );

  ThrowErrorMsgIf(duplicate && existing_rank != entity_rank,
    "For args: cell_topolgy " << cell_topology.getName() << " and entity_rank " << entity_rank << ", " <<
    "previously declared rank = " << existing_rank );

  if (! duplicate) {
    std::string part_name = std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + std::string(cell_topology.getName());

    ThrowErrorMsgIf(get_part(part_name) != 0, "Cannot register topology with same name as existing part '" << cell_topology.getName() << "'" );

    Part &part = declare_internal_part(part_name, entity_rank);
    m_cellTopologyPartEntityRankMap[cell_topology] = CellTopologyPartEntityRankMap::mapped_type(&part, entity_rank);

    assign_cell_topology(m_partCellTopologyVector, part.mesh_meta_data_ordinal(), cell_topology);
  }
  //check_topo_db();
}


CellTopology
MetaData::get_cell_topology(
  const std::string &   topology_name) const
{
  std::string part_name = convert_to_internal_name(std::string("FEM_ROOT_CELL_TOPOLOGY_PART_") + topology_name);

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

EntityRank MetaData::get_entity_rank(
  const CellTopology       cell_topology) const
{
  CellTopologyPartEntityRankMap::const_iterator it = m_cellTopologyPartEntityRankMap.find(cell_topology);
  if (it == m_cellTopologyPartEntityRankMap.end())
    return INVALID_RANK;
  else
    return (*it).second.second;
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
    const PartVector & intersect = p.intersection_of();

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
    const unsigned intersect_size = static_cast<unsigned>(intersect.size());
    b.pack<unsigned>( intersect_size );
    for ( j = intersect.begin() ; j != intersect.end() ; ++j ) {
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
    const PartVector & intersect = p.intersection_of();
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

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == intersect.size();
    }
    for ( j = intersect.begin() ; ok && j != intersect.end() ; ++j ) {
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

bool is_cell_topology_root_part(const Part & part) {
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
void set_cell_topology(
  Part &                        part,
  CellTopology             cell_topology)
{
  MetaData& meta = MetaData::get(part);

  ThrowRequireMsg(meta.is_initialized(),"set_cell_topology: initialize() must be called before this function");

  Part &root_part = meta.get_cell_topology_root_part(cell_topology);
  meta.declare_part_subset(root_part, part);
}

std::vector<std::string>
entity_rank_names( size_t spatial_dimension )
{
  ThrowInvalidArgMsgIf( spatial_dimension < 1 || 3 < spatial_dimension,
                        "Invalid spatial dimension = " << spatial_dimension );

  std::vector< std::string > names ;

  names.reserve( spatial_dimension + 1 );

  names.push_back( std::string( "NODE" ) );

  if ( 1 < spatial_dimension ) { names.push_back( std::string("EDGE") ); }
  if ( 2 < spatial_dimension ) { names.push_back( std::string("FACE") ); }

  names.push_back( std::string("ELEMENT") );

  return names ;
}


CellTopology
get_cell_topology(
  const Bucket &                bucket)
{
  const BulkData   & bulk_data = BulkData::get(bucket);
  const MetaData   & meta_data = MetaData::get(bulk_data);
  const PartVector & all_parts = meta_data.get_parts();

  CellTopology cell_topology;

  const std::pair< const unsigned *, const unsigned * > supersets = bucket.superset_part_ordinals();

  if (supersets.first != supersets.second) {
    const Part *first_found_part = 0;

    for ( const unsigned * it = supersets.first ; it != supersets.second ; ++it ) {

      const Part & part = * all_parts[*it] ;

      if ( part.primary_entity_rank() == bucket.entity_rank() ) {

        CellTopology top = meta_data.get_cell_topology( part );

        if ( ! cell_topology.getCellTopologyData() ) {
          cell_topology = top ;

          if (!first_found_part)
            first_found_part = &part;
        }
        else {
          ThrowErrorMsgIf( top.getCellTopologyData() && top != cell_topology,
            "Cell topology is ambiguously defined. It is defined as " << cell_topology.getName() <<
            " on part " << first_found_part->name() << " and as " << top.getName() << " on its superset part " << part.name() );
        }
      }
    }
  }

  return cell_topology ;
}

CellTopology get_cell_topology(const Entity &entity)
{
  return get_cell_topology(entity.bucket());
}

} // namespace mesh
} // namespace stk

