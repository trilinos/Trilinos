/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_BucketImpl_hpp
#define stk_mesh_BucketImpl_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

class BucketImpl {
  public:
  ~BucketImpl();

  struct DataMap {
    typedef FieldBase::Restriction::size_type size_type ;
    const size_type * m_stride ;
    size_type         m_base ;
    size_type         m_size ;
  };

  BucketImpl( BulkData & arg_mesh ,
          unsigned          arg_entity_rank ,
          const unsigned  * arg_key ,
          size_t            arg_alloc_size ,
          size_t            arg_capacity ,
          DataMap         * arg_field_map ,
          Entity         ** arg_entity_array );

  //
  // External interface:
  //
  BulkData & mesh() const { return m_mesh ; }
  unsigned entity_rank() const { return m_entity_rank ; }
  const unsigned * key() const { return m_key ; }
  std::pair<const unsigned *, const unsigned *>
    superset_part_ordinals() const
    {
      return std::pair<const unsigned *, const unsigned *>
             ( m_key + 1 , m_key + m_key[0] );
    }
  unsigned allocation_size() const { return m_alloc_size ; }
  size_t capacity() const { return m_capacity ; }
  size_t size() const { return m_size ; }
  Entity & operator[] ( size_t i ) const { return *(m_entities[i]) ; }
  unsigned field_data_size(const FieldBase & field) const
  {
    return m_field_map[ field.mesh_meta_data_ordinal() ].m_size;
  }
  const FieldBase::Restriction::size_type * field_data_stride( const FieldBase & field ) const
  {
    return m_field_map[ field.mesh_meta_data_ordinal() ].m_stride;
  }
  unsigned char * field_data_location( const FieldBase & field, const Entity & entity ) const
  {
    return field_data_location_impl( field.mesh_meta_data_ordinal(), entity.bucket_ordinal() );
  }
  unsigned char * field_data_location( const FieldBase & field ) const
  {
    unsigned int zero_ordinal = 0;
    return field_data_location_impl( field.mesh_meta_data_ordinal(), zero_ordinal );
  }

  //
  // Internal interface:
  //
  void increment_size() { ++m_size ; }
  void decrement_size() { --m_size ; }
  void replace_entity(unsigned entity_ordinal, Entity * entity ) { m_entities[entity_ordinal] = entity ; }
  void update_state();

  template< class field_type >
  typename FieldTraits< field_type >::data_type *
  field_data( const field_type & f , const unsigned & entity_ordinal ) const
  {
    typedef typename FieldTraits< field_type >::data_type * data_p ;
    return (data_p)(field_data_location_impl(f.mesh_meta_data_ordinal(),entity_ordinal));
  }

  // BucketKey key = ( part-count , { part-ordinals } , counter )
  //  key[ key[0] ] == counter
  unsigned bucket_counter() const { return m_key[ *m_key ]; }

  Bucket * last_bucket_in_family();
  Bucket * first_bucket_in_family();
  void set_last_bucket_in_family( Bucket * last_bucket );
  void set_first_bucket_in_family( Bucket * first_bucket );
  DataMap * get_field_map();
  void zero_fields( unsigned i_dst );
  void replace_fields( unsigned i_dst , Bucket & k_src , unsigned i_src );
  void set_bucket_family_pointer( Bucket * bucket ) { m_bucket = bucket; }
  const Bucket * get_bucket_family_pointer() const { return m_bucket; }
  bool equivalent( const BucketImpl& other_bucket ) const {
    return m_bucket == other_bucket.m_bucket;
  }

  private:
  BucketImpl();

  BulkData             & m_mesh ;        // Where this bucket resides
  const unsigned         m_entity_rank ; // Type of entities for this bucket
  const unsigned * const m_key ;         // Unique key in the bulk data
  const size_t           m_alloc_size ;  // Allocation size of this bucket
  const size_t           m_capacity ;    // Capacity for entities
  size_t                 m_size ;        // Number of entities
  Bucket               * m_bucket ;      // Pointer to head of bucket family, but head points to tail
  DataMap        * const m_field_map ;   // Field value data map, shared
  Entity        ** const m_entities ;    // Array of entity pointers,
                                         // beginning of field value memory.

  unsigned char * field_data_location_impl( const unsigned & field_ordinal, const unsigned & entity_ordinal ) const
  {
    typedef unsigned char * byte_p ;
    const DataMap & data_map = m_field_map[ field_ordinal ];
    unsigned char * ptr = NULL;
    if ( data_map.m_size ) {
      ptr = ((byte_p)(m_entities) + data_map.m_base + data_map.m_size * entity_ordinal );
    }
    return ptr ;
  }
  Bucket * last_bucket_in_family_impl();
};



} // namespace impl
} // namespace mesh
} // namespace stk


#endif // stk_mesh_BucketImpl_hpp

