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

#include <stk_mesh/base/BulkData.hpp>


#include <stk_mesh/baseImpl/FieldRepository.hpp>

namespace stk {
namespace mesh {

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

MetaData::MetaData(const std::vector<std::string>& entity_rank_names)
  : m_commit( false ),
    m_part_repo( this ),
    m_attributes(),
    m_universal_part( NULL ),
    m_owns_part( NULL ),
    m_shares_part( NULL ),
    m_field_repo(),
    m_field_relations( ),
    m_properties( ),
    m_entity_rank_names( entity_rank_names )
{
  ThrowErrorMsgIf( entity_rank_names.empty(), "entity ranks empty" );

  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_part( std::string("{OWNS}") );
  m_shares_part = & declare_part( std::string("{SHARES}") );
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
    m_entity_rank_names( )
{
  // Declare the predefined parts

  m_universal_part = m_part_repo.universal_part();
  m_owns_part = & declare_part( std::string("{OWNS}") );
  m_shares_part = & declare_part( std::string("{SHARES}") );
}

//----------------------------------------------------------------------

void MetaData::set_entity_rank_names(const std::vector<std::string> &entity_rank_names)
{
  ThrowErrorMsgIf( entity_rank_names.empty(), "entity ranks empty" );

  m_entity_rank_names = entity_rank_names;
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


Part & MetaData::declare_part( const std::string & p_name , EntityRank rank )
{
  require_not_committed();
  require_valid_entity_rank(rank);

  return *m_part_repo.declare_part( p_name , rank );
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
  static const char method[] = "stk::mesh::MetaData::declare_part_subset" ;

  require_not_committed();
  require_same_mesh_meta_data( MetaData::get(superset) );
  require_same_mesh_meta_data( MetaData::get(subset) );
  require_not_relation_target( &superset );
  require_not_relation_target( &subset );

  m_part_repo.declare_subset( superset, subset );

  // The new superset / subset relationship can cause a
  // field restriction to become incompatible or redundant.
  m_field_repo.verify_and_clean_restrictions(method, m_part_repo.get_all_parts());
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
  const unsigned * arg_stride )
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
      arg_stride
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

} // namespace mesh
} // namespace stk

