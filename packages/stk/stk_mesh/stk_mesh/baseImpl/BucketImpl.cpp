/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <stk_mesh/baseImpl/BucketImpl.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------

namespace {

void memory_copy( unsigned char * dst , const unsigned char * src , unsigned n )
{ std::memcpy( dst , src , n ); }


void memory_zero( unsigned char * dst , unsigned n )
{ std::memset( dst , 0 , n ); }

} // namespace

//----------------------------------------------------------------------

void BucketImpl::update_state()
{
  const MetaData & meta = MetaData::get(m_mesh);
  const std::vector<FieldBase*> & field_set = meta.get_fields();

  for ( unsigned i = 0 ; i < field_set.size() ; ) {

    DataMap * const tmp = &m_field_map[0] + i ;
    const FieldBase & field = * field_set[i] ;
    const unsigned num_state = field.number_of_states();
    i += num_state ;

    if ( 1 < num_state && tmp->m_size ) {
      unsigned offset[ MaximumFieldStates ] ;

      for ( unsigned j = 0 ; j < num_state ; ++j ) {
        offset[j] = tmp[j].m_base ;
      }

      for ( unsigned j = 0 ; j < num_state ; ++j ) {
        const unsigned j_new = ( j + num_state - 1 ) % num_state ;
        tmp[j_new].m_base = offset[j] ;
      }
    }
  }
}

//----------------------------------------------------------------------
// Every bucket in the family points to the first bucket,
// except the first bucket which points to the last bucket.

Bucket * BucketImpl::last_bucket_in_family() const
{
  Bucket * last = last_bucket_in_family_impl();

  ThrowRequireMsg( NULL != last, "Last is NULL");
  ThrowRequireMsg( last->size() != 0, "Last bucket is empty");

  return last ;
}

Bucket * BucketImpl::last_bucket_in_family_impl() const
{
  bool this_is_first_bucket_in_family = (bucket_counter() == 0);

  Bucket * last = NULL;

  if (this_is_first_bucket_in_family) {
    last = m_bucket;
  } else {
    last = m_bucket->m_bucketImpl.m_bucket;
  }

  return last;
}

//----------------------------------------------------------------------

Bucket * BucketImpl::first_bucket_in_family() const
{
  return last_bucket_in_family_impl()->m_bucketImpl.m_bucket;
}

//----------------------------------------------------------------------

void BucketImpl::set_last_bucket_in_family( Bucket * last_bucket )
{
  Bucket * last = last_bucket_in_family_impl();
  Bucket * first = last->m_bucketImpl.m_bucket;
  first->m_bucketImpl.m_bucket = last_bucket;
}

//----------------------------------------------------------------------

void BucketImpl::set_first_bucket_in_family( Bucket * first_bucket )
{
  m_bucket = first_bucket;
}

//----------------------------------------------------------------------

BucketImpl::DataMap * BucketImpl::get_field_map()
{
  return &m_field_map[0];
}

//----------------------------------------------------------------------

void BucketImpl::initialize_fields( unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const p = &m_field_data[0];
  const DataMap *       i = &m_field_map[0];
  const DataMap * const e = i + field_set.size();

  for (std::vector<FieldBase*>::const_iterator field_iter=field_set.begin() ;
       i != e ; ++i, ++field_iter ) {

    if (i->m_size == 0) continue;

    const unsigned char* init_val = reinterpret_cast<const unsigned char*>((*field_iter)->get_initial_value());
    if (init_val != NULL) {
      memory_copy( p + i->m_base + i->m_size * i_dst , init_val, i->m_size );
    }
    else {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void BucketImpl::replace_fields( unsigned i_dst , Bucket & k_src , unsigned i_src )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const s = &(k_src.m_bucketImpl.m_field_data[0]);
  unsigned char * const d = &m_field_data[0];
  const DataMap *       j = &(k_src.m_bucketImpl.m_field_map[0]);
  const DataMap *       i = &m_field_map[0];
  const DataMap * const e = i + field_set.size();

  for (std::vector<FieldBase*>::const_iterator field_iter=field_set.begin() ;
       i != e ; ++i , ++j, ++field_iter ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        ThrowErrorMsgIf( i->m_size != j->m_size,
            "Incompatible field sizes: " << i->m_size << " != " << j->m_size );

        memory_copy( d + i->m_base + i->m_size * i_dst ,
                     s + j->m_base + j->m_size * i_src , i->m_size );
      }
      else {
        const unsigned char* init_val = reinterpret_cast<const unsigned char*>((*field_iter)->get_initial_value());
        if (init_val != NULL) {
          memory_copy( d + i->m_base + i->m_size * i_dst ,
                       init_val, i->m_size );
        }
        else {
          memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
        }
      }
    }
  }
}

//----------------------------------------------------------------------
namespace {
inline unsigned align( size_t nb )
{
  enum { BYTE_ALIGN = 16 };
  const unsigned gap = nb % BYTE_ALIGN ;
  if ( gap ) { nb += BYTE_ALIGN - gap ; }
  return nb ;
}

const FieldBase::Restriction & empty_field_restriction()
{
  static const FieldBase::Restriction empty ;
  return empty ;
}

const FieldBase::Restriction & find_restriction( const FieldBase & field ,
                                          EntityRank erank ,
                                          const unsigned num_part_ord ,
                                          const unsigned part_ord[] )
{
  const FieldBase::Restriction & empty = empty_field_restriction();
  const FieldBase::Restriction * restriction = & empty ;

  const std::vector<FieldBase::Restriction> & restr_vec = field.restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator iend = restr_vec.end();
        std::vector<FieldBase::Restriction>::const_iterator ibeg = restr_vec.begin();

  for ( PartOrdinal i = 0 ; i < num_part_ord && iend != ibeg ; ++i ) {

    const FieldRestriction restr(erank,part_ord[i]);

    //lower_bound returns an iterator to either the insertion point for the
    //'restr' argument, or to a matching restriction.
    //It only returns the 'end' iterator if 'restr' is past the end of the
    //vector of restrictions being searched.
    //This depends on the input part ordinals being sorted, and on the restriction
    //vector being sorted by part ordinal.

    ibeg = std::lower_bound( ibeg , iend , restr );

    if ( (iend != ibeg) && (*ibeg == restr) ) {
      if ( restriction == & empty ) { restriction = & *ibeg ; }

      if ( ibeg->not_equal_stride(*restriction) ) {

        Part & p_old = MetaData::get(field).get_part( ibeg->part_ordinal() );
        Part & p_new = MetaData::get(field).get_part( restriction->part_ordinal() );

        std::ostringstream msg ;
        msg << " FAILED WITH INCOMPATIBLE DIMENSIONS FOR " ;
        msg << field ;
        msg << " Part[" << p_old.name() ;
        msg << "] and Part[" << p_new.name() ;
        msg << "]" ;

        ThrowErrorMsg( msg.str() );
      }
    }
  }

  const std::vector<FieldBase::Restriction> & sel_res = field.selector_restrictions();
  std::pair<const unsigned*,const unsigned*> bucket_part_range = std::make_pair(part_ord, part_ord+num_part_ord);
  for(std::vector<FieldBase::Restriction>::const_iterator it=sel_res.begin(), it_end=sel_res.end(); it != it_end; ++it) {
    const Selector& selector = it->selector();
    if (it->entity_rank() == erank && selector.apply(bucket_part_range)) {
      if (restriction == &empty) {
        restriction = &*it;
      }
      if (it->not_equal_stride(*restriction)) {
        ThrowErrorMsg("find_restriction calculation failed with different field-restriction selectors giving incompatible sizes.");
      }
    }
  }

  return *restriction ;
}
} // namespace
//----------------------------------------------------------------------

BucketImpl::BucketImpl( BulkData & arg_mesh,
                        EntityRank arg_entity_rank,
                        const std::vector<unsigned> & arg_key,
                        size_t arg_capacity
                      )
  : m_mesh(arg_mesh)
  , m_entity_rank(arg_entity_rank)
  , m_key(arg_key)
  , m_capacity(arg_capacity)
  , m_size(0)
  , m_bucket(NULL)
  , m_field_map( m_mesh.mesh_meta_data().get_fields().size()+1)
  , m_entities(arg_capacity)
  , m_field_data(0)
{

  //calculate the size of the field_data

  const std::vector< FieldBase * > & field_set =
    arg_mesh.mesh_meta_data().get_fields();

  const size_t num_fields = field_set.size();

  size_t field_data_size = 0;

  if (arg_capacity != 0) {
    for ( size_t i = 0; i<num_fields; ++i) {
      const FieldBase  & field = * field_set[i] ;
      unsigned num_bytes_per_entity = 0 ;

      const FieldBase::Restriction & restriction =
        find_restriction( field, arg_entity_rank, m_key[0]-1, &m_key[1]);

      if ( restriction.dimension() > 0 ) { // Exists

        const unsigned type_stride = field.data_traits().stride_of ;
        const unsigned field_rank  = field.rank();

        num_bytes_per_entity = type_stride *
          ( field_rank ? restriction.stride( field_rank - 1 ) : 1 );
      }
      m_field_map[i].m_base = field_data_size ;
      m_field_map[i].m_size = num_bytes_per_entity ;
      m_field_map[i].m_stride = &restriction.stride(0);

      field_data_size += align( num_bytes_per_entity * m_capacity );
    }
    m_field_map[ num_fields ].m_base  = field_data_size ;
    m_field_map[ num_fields ].m_size = 0 ;
    m_field_map[ num_fields ].m_stride = NULL ;
  }
  else { //nil bucket

    FieldBase::Restriction::size_type empty_stride[ MaximumFieldDimension ];
    Copy<MaximumFieldDimension>( empty_stride , FieldBase::Restriction::size_type(0) );

    for ( size_t i = 0; i<num_fields; ++i) {
      m_field_map[i].m_base = 0 ;
      m_field_map[i].m_size = 0 ;
      m_field_map[i].m_stride = empty_stride;
    }
    m_field_map[ num_fields ].m_base   = 0 ;
    m_field_map[ num_fields ].m_size   = 0 ;
    m_field_map[ num_fields ].m_stride = NULL ;
  }

  //allocate space for the fields
  m_field_data.resize(field_data_size);

}

} // namespace impl
} // namespace mesh
} // namespace stk
