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

#include <boost/pool/pool_alloc.hpp>
//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

class BucketImpl {
  public:

  struct DataMap {
    typedef FieldBase::Restriction::size_type size_type ;
    const size_type * m_stride ;
    size_type         m_base ;
    size_type         m_size ;
  };

  BucketImpl( BulkData & arg_mesh ,
              EntityRank arg_entity_rank,
              const std::vector<unsigned> & arg_key,
              size_t arg_capacity
            );

  //
  // External interface:
  //
  BulkData & mesh() const { return m_mesh ; }
  unsigned entity_rank() const { return m_entity_rank ; }
  const unsigned * key() const { return &m_key[0] ; }
  const std::vector<unsigned> & key_vector() const { return m_key; }

  std::pair<const unsigned *, const unsigned *>
    superset_part_ordinals() const
    {
      return std::pair<const unsigned *, const unsigned *>
             ( key() + 1 , key() + key()[0] );
    }
  unsigned allocation_size() const { return 0 ; }
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

  /** Experimental method that skips the if-test for field existence on the bucket.
    This method should only be called if the caller is sure that the field exists on the bucket.
  */
  unsigned char * fast_field_data_location( const FieldBase & field, unsigned ordinal ) const
  {
    return fast_field_data_location_impl( field.mesh_meta_data_ordinal(), ordinal );
  }
  unsigned char * field_data_location( const FieldBase & field, unsigned ordinal ) const
  {
    return field_data_location_impl( field.mesh_meta_data_ordinal(), ordinal );
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
    return reinterpret_cast<data_p>(field_data_location_impl(f.mesh_meta_data_ordinal(),entity_ordinal));
  }

  // BucketKey key = ( part-count , { part-ordinals } , counter )
  //  key[ key[0] ] == counter
  unsigned bucket_counter() const { return m_key[ m_key[0] ]; }

  Bucket * last_bucket_in_family() const;
  Bucket * first_bucket_in_family() const;
  void set_last_bucket_in_family( Bucket * last_bucket );
  void set_first_bucket_in_family( Bucket * first_bucket );
  DataMap * get_field_map();
  void initialize_fields( unsigned i_dst );
  void replace_fields( unsigned i_dst , Bucket & k_src , unsigned i_src );
  void set_bucket_family_pointer( Bucket * bucket ) { m_bucket = bucket; }
  const Bucket * get_bucket_family_pointer() const { return m_bucket; }

  bool equivalent( const BucketImpl& other_bucket ) const {
    return first_bucket_in_family() == other_bucket.first_bucket_in_family();
  }

  Entity*const* begin() const { return &m_entities[0]; }
  Entity*const* end() const { return &m_entities[0] + m_size; }

  ~BucketImpl() { delete [] m_field_data; }

  private:
  BucketImpl();

  BulkData             & m_mesh ;        // Where this bucket resides
  const EntityRank       m_entity_rank ; // Type of entities for this bucket
  std::vector<unsigned>  m_key ;
  const size_t           m_capacity ;    // Capacity for entities
  size_t                 m_size ;        // Number of entities
  Bucket               * m_bucket ;      // Pointer to head of bucket family, but head points to tail
  std::vector<DataMap>   m_field_map ;   // Field value data map, shared
  std::vector<Entity*>   m_entities ;    // Array of entity pointers,
                                         // beginning of field value memory.
  unsigned char* m_field_data;
  unsigned char* m_field_data_end;

  unsigned char * field_data_location_impl( const unsigned & field_ordinal, const unsigned & entity_ordinal ) const
  {
    typedef unsigned char * byte_p ;
    const DataMap & data_map = m_field_map[ field_ordinal ];
    unsigned char * ptr = NULL;
    if ( data_map.m_size ) {
      ptr = const_cast<unsigned char*>(m_field_data) + data_map.m_base + data_map.m_size * entity_ordinal;
      ThrowAssert(ptr < m_field_data_end);
    }
    return ptr ;
  }
  unsigned char * fast_field_data_location_impl( const unsigned & field_ordinal, const unsigned & entity_ordinal ) const
  {
    typedef unsigned char * byte_p ;
    const DataMap & data_map = m_field_map[ field_ordinal ];
    ThrowAssertMsg(data_map.m_size>0,"Field doesn't exist on bucket.");
    return const_cast<unsigned char*>(m_field_data) + data_map.m_base + data_map.m_size * entity_ordinal;
  }
  Bucket * last_bucket_in_family_impl() const;
};



} // namespace impl
} // namespace mesh
} // namespace stk


#endif // stk_mesh_BucketImpl_hpp

