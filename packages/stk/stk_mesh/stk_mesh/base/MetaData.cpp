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

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void MetaData::assert_not_committed( const char * method ) const
{
  if ( m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh MetaData has been committed." );
    throw std::logic_error( msg );
  }
}

void MetaData::assert_committed( const char * method ) const
{
  if ( ! m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh MetaData has not been committed." );
    throw std::logic_error( msg );
  }
}

void MetaData::assert_same_mesh_meta_data( const char * method ,
                                           const MetaData & rhs ) const
{
  if ( this != & rhs ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: Different mesh_meta_data." );
    throw std::logic_error( msg );
  }
}

void MetaData::assert_entity_type( const char * method ,
                                   unsigned rank ) const
{
  if ( m_entity_type_names.size() <= rank ) {
    std::ostringstream msg ;
    msg << method ;
    msg << " FAILED: entity_type( " << rank ;
    msg << " ) >= maximum_value( " << m_entity_type_names.size();
    msg << " )" ;
    throw std::logic_error( msg.str() );
  }
}

//----------------------------------------------------------------------

MetaData::MetaData(const std::vector<std::string>& entity_type_names)
  : m_commit( false ),
    m_universal_part( NULL ),
    m_uses_part( NULL ),
    m_owns_part( NULL ),
    m_fields( ),
    m_field_relations( ),
    m_properties( ),
    m_entity_type_names( entity_type_names )
{
  if ( entity_type_names.empty() ) {
    std::string msg( "stk::mesh::MetaData constructor FAILED: no entity types" );
    throw std::runtime_error( msg );
  }

  // Declare the predefined parts

  m_universal_part = new Part( this );
  m_uses_part = & declare_part( std::string("{USES}") );
  m_owns_part = & declare_part( std::string("{OWNS}") );

  declare_part_subset( * m_uses_part , * m_owns_part );
}

//----------------------------------------------------------------------

const std::string& MetaData::entity_type_name( unsigned ent_type ) const
{
  if (ent_type >= m_entity_type_names.size()) {
    std::ostringstream msg;
    msg << "Error in MetaData::entity_type_name: entity-type (" << ent_type
        << ") out of range. Must be in range [0 .. " << m_entity_type_names.size()
        << ").";
    throw std::runtime_error( msg.str() );
  }

  return m_entity_type_names[ent_type];
}

//----------------------------------------------------------------------

Part * MetaData::get_part( const std::string & p_name ,
                           const char * required_by ) const
{
  const PartVector & all_parts = m_universal_part->subsets();

  Part * const p = find( all_parts , p_name );

  if ( required_by && NULL == p ) { // ERROR
    static const char method[] = "stk::mesh::BulkData::get_part" ;
    std::string msg ;
    msg.append( method )
       .append( "( " )
       .append( p_name )
       .append( " , " )
       .append( required_by )
       .append( " ) FAILED to find part" );
    throw std::runtime_error( msg );
  }

  return p ;
}

Part & MetaData::declare_part( const std::string & p_name )
{
  static const char method[] = "stk::mesh::MetaData::declare_part" ;

  const unsigned rank = std::numeric_limits<unsigned>::max();

  assert_not_committed( method );

  return m_universal_part->declare_part( p_name , rank );
}

Part & MetaData::declare_part( const std::string & p_name , EntityType rank )
{
  static const char method[] = "stk::mesh::MetaData::declare_part" ;

  assert_not_committed( method );
  assert_entity_type( method , rank );

  return m_universal_part->declare_part( p_name , rank );
}

namespace {

void assert_not_relation_target(
  const char * const method ,
  const Part * const part )
{
  std::vector<PartRelation>::const_iterator i_end = part->relations().end();
  std::vector<PartRelation>::const_iterator i     = part->relations().begin();
  for ( ; i != i_end ; ++i ) {
    if ( part == i->m_target ) {
      std::string msg ;
      msg.append( method );
      msg.append( "(...) FAILED Requirement that Part[" );
      msg.append( part->name() );
      msg.append( "] is not a PartRelation target" );
      throw std::runtime_error(msg);
    }
  }
}

}

Part & MetaData::declare_part( const PartVector & part_intersect )
{
  static const char method[] = "stk::mesh::MetaData::declare_part" ;

  assert_not_committed( method );

  for ( PartVector::const_iterator
        i = part_intersect.begin() ; i != part_intersect.end() ; ++i ) {
    assert_not_relation_target( method , *i );
  }

  return m_universal_part->declare_part( part_intersect );
}

void MetaData::declare_part_subset( Part & superset , Part & subset )
{
  static const char method[] = "stk::mesh::MetaData::declare_part_subset" ;

  assert_not_committed( method );
  assert_same_mesh_meta_data( method , superset.mesh_meta_data() );
  assert_same_mesh_meta_data( method , subset.mesh_meta_data() );
  assert_not_relation_target( method , & superset );
  assert_not_relation_target( method , & subset );

  superset.declare_subset( subset );

  // The new superset / subset relationship can cause a
  // field restriction to become incompatible or redundant.

  const std::vector<Part*> & all_parts = m_universal_part->subsets();

  for ( std::vector<FieldBase *>::iterator
        f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    (*f)->verify_and_clean_restrictions( method , all_parts );
  }
}

void MetaData::declare_part_relation(
  Part & root_part ,
  relation_stencil_ptr stencil ,
  Part & target_part )
{
  static const char method[] = "stk::mesh::MetaData::declare_part_relation" ;

  assert_not_committed( method );
  assert_not_relation_target( method , & root_part );

  if (!stencil) {
    std::string msg ;
    msg.append( method );
    msg.append( "stencil function pointer cannott be NULL" );
    throw std::runtime_error( msg );
  }

  if ( 0 != target_part.subsets().size() ||
       0 != target_part.intersection_of().size() ||
       1 != target_part.supersets().size() ) {
    std::string msg ;
    msg.append( method );
    msg.append( ": FAILED Requirement that target Part[" );
    msg.append( target_part.name() );
    msg.append( "] is not a superset or subset" );
    throw std::runtime_error( msg );
  }

  PartRelation tmp ;
  tmp.m_root = & root_part ;
  tmp.m_target = & target_part ;
  tmp.m_function = stencil ;

  root_part.m_relations.push_back( tmp );
  target_part.m_relations.push_back( tmp );
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
  static const char method[] = "std::mesh::MetaData::declare_field" ;

  assert_not_committed( method );

  return FieldBase::declare_field( arg_name , arg_traits ,
                                   arg_rank , arg_dim_tags ,
                                   arg_num_states ,
                                   this , m_fields );
}

void MetaData::declare_field_restriction(
  FieldBase      & arg_field ,
  unsigned         arg_entity_type ,
  const Part     & arg_part ,
  const unsigned * arg_stride )
{
  static const char method[] =
    "std::mesh::MetaData::declare_field_restriction" ;

  assert_not_committed( method );
  assert_same_mesh_meta_data( method , arg_field.mesh_meta_data() );
  assert_same_mesh_meta_data( method , arg_part.mesh_meta_data() );

  arg_field.insert_restriction( method, arg_entity_type, arg_part, arg_stride);

  arg_field.verify_and_clean_restrictions( method, get_parts() );
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

void MetaData::declare_field_lock_relation(
  FieldBase & pointer_field ,
  relation_stencil_ptr stencil )
{
  FieldRelation tmp ;
  tmp.m_root   = & pointer_field ;
  tmp.m_target = 0 ;
  tmp.m_function = stencil ;

  m_field_relations.push_back( tmp );
}

//----------------------------------------------------------------------

void MetaData::commit()
{
  static const char method[] = "stk::mesh::MetaData::commit" ;

  assert_not_committed( method );

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

  // Destroy the fields, used 'new' to allocate so now use 'delete'

  try {
    std::vector<FieldBase * >::iterator j = m_fields.begin();

    for ( ; j != m_fields.end() ; ++j ) { delete *j ; }

    m_fields.clear();
  } catch(...) {}

  try {
    delete m_universal_part ;
  } catch(...){}
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
  static const char method[] = "stk::mesh::verify_parallel_consistency(MetaData)" ;

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

  if ( ! ok[0] || ! ok[1] ) {
    std::ostringstream msg ;
    msg << "P" << p_rank ;
    msg << ": " << method ;
    msg << " : FAILED for:" ;
    if ( ! ok[0] ) { msg << " Parts" ; }
    if ( ! ok[1] ) { msg << " Fields" ; }
    throw std::logic_error( msg.str() );
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

