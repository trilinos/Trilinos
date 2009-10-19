
#ifndef stk_mesh_FieldData_hpp
#define stk_mesh_FieldData_hpp

//----------------------------------------------------------------------

#include <Shards_Array.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_field_data
 *  \{
 */

/** \class BucketArray
 *  \brief  \ref stk::mesh::Field "Field" data \ref shards::Array "Array"
 *          for a given array field and bucket
 */
template< class FieldType > struct BucketArray {};

/** \class EntityArray
 *  \brief  \ref stk::mesh::Field "Field" data \ref shards::Array "Array"
 *          for a given array field and entity
 */
template< class FieldType > struct EntityArray {};

//----------------------------------------------------------------------
/** \brief  Check for existence of field data.
 *
 *  \exception std::runtime_error
 *     Thrown if required_by != NULL and the field data does not exist.
 */
bool field_data_valid( const FieldBase & f ,
                       const Bucket & k ,
                       unsigned ord = 0,
                       const char * required_by = NULL );

/** \brief  Check for existence of field data.
 *
 *  \exception std::runtime_error
 *     Thrown if required_by != NULL and the field data does not exist.
 */
inline
bool field_data_valid( const FieldBase & f ,
                       const Entity & e ,
                       const char * required_by = NULL )
{ return field_data_valid( f, e.bucket(), e.bucket_ordinal(), required_by ); }

//----------------------------------------------------------------------

/** \brief  Size, in bytes, of the field data for each entity */
inline
unsigned field_data_size( const FieldBase & f , const Bucket & k )
{
  const Bucket::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];
  return pd.m_size ;
}

/** \brief  Size, in bytes, of the field data for each entity */
inline
unsigned field_data_size( const FieldBase & f , const Entity & e )
{ return field_data_size( f , e.bucket() ); }

//----------------------------------------------------------------------


/** \brief  Pointer to the field data array */
template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Bucket::iterator &i ) 
{
  typedef unsigned char * byte_p ;
  typedef typename FieldTraits< field_type >::data_type * data_p ;

  data_p ptr = NULL ;

  {
    const Bucket          & b  = * i.m_bucket_ptr ;
    const Bucket::DataMap & pd = b.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      ptr = (data_p)( ((byte_p)(b.m_entities)) +
                      pd.m_base + pd.m_size * i.m_current_entity );
    }
  }
  return ptr ;

}


/** \brief  Pointer to the field data array */
template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Entity & e )
{
  typedef unsigned char * byte_p ;
  typedef typename FieldTraits< field_type >::data_type * data_p ;

  data_p ptr = NULL ;

  {
    const Bucket & b = e.bucket();

    const Bucket::DataMap & pd = b.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      ptr = (data_p)( ((byte_p)(b.m_entities)) +
                      pd.m_base + pd.m_size * e.bucket_ordinal() );
    }
  }
  return ptr ;
}

//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

void throw_field_data_array( const FieldBase & , unsigned );

template< typename ScalarType >
struct EntityArray< Field<ScalarType,void,void,void,void,void,void,void> >
  : public shards::Array<ScalarType,shards::RankZero,void,void,void,void,void,void,void>
{
  typedef Field<ScalarType,void,void,void,void,void,void,void> field_type ;
  typedef shards::Array<ScalarType,shards::RankZero,void,void,void,void,void,void,void>
  array_type ;

  EntityArray( const field_type & f , const Entity & e )
    : array_type( field_data( f , e ) ) {}

private:
  EntityArray();
  EntityArray( const EntityArray & );
  EntityArray & operator = ( const EntityArray & );
};
#endif /* DOXYGEN_COMPILE */

/** \brief  \ref stk::mesh::Field "Field" data \ref shards::Array "Array"
 *          for a given array field and entity
 */
template< typename ScalarType ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct EntityArray< Field<ScalarType,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
  : public shards::Array<ScalarType,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
{
#ifndef DOXYGEN_COMPILE
private:
  typedef unsigned char * byte_p ;
  EntityArray();
  EntityArray( const EntityArray & );
  EntityArray & operator = ( const EntityArray & );
public:

  typedef Field<ScalarType,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;

  typedef
  shards::Array<ScalarType,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
  array_type ;

  EntityArray( const field_type & f , const Entity & e ) : array_type()
  {
    const Bucket          & b = e.bucket();
    const Bucket::DataMap & pd = b.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      array_type::assign_stride( 
        (ScalarType*)( (byte_p)(b.m_entities) +
                       pd.m_base + pd.m_size * e.bucket_ordinal() ),
        pd.m_stride );
    }
  }
#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------
/** \brief Implement ArrayDimTag for the entity count dimension
 *         of a BucketArray.
 * 
 *  \ref stk_mesh_field_data "array field data" obtained from a
 *  \ref stk::mesh::Bucket "bucket".
 */ 
struct EntityDimension : public shards::ArrayDimTag {
   
  const char * name() const ;
   
  static const EntityDimension & tag(); ///< Singleton
 
private: 
  EntityDimension() {} 
  EntityDimension( const EntityDimension & ); 
  EntityDimension & operator = ( const EntityDimension & );
};


/** \brief  \ref stk::mesh::Field "Field" data \ref shards::Array "Array"
 *          for a given scalar field and bucket
 */
template< typename ScalarType >
struct BucketArray< Field<ScalarType,void,void,void,void,void,void,void> >
  : public
shards::Array<ScalarType,shards::FortranOrder,EntityDimension,void,void,void,void,void,void> 
{
#ifndef DOXYGEN_COMPILE
private:
  typedef unsigned char * byte_p ;
  BucketArray();
  BucketArray( const BucketArray & );
  BucketArray & operator = ( const BucketArray & );

public:

  typedef Field<ScalarType,void,void,void,void,void,void,void> field_type ;
  typedef
  shards::Array<ScalarType,shards::FortranOrder,EntityDimension,void,void,void,void,void,void> 
  array_type ;

  BucketArray( const field_type & f , const Bucket & k )
  {
    const Bucket::DataMap & pd = k.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      array_type::assign( (ScalarType*)( (byte_p)(k.m_entities) + pd.m_base ) ,
                          k.size() );
    }
  }

  BucketArray( const field_type & f,
               const Bucket::iterator & i,
               const Bucket::iterator & j) 
  {
    const Bucket::DataMap & pd =
      i.m_bucket_ptr->m_field_map[ f.mesh_meta_data_ordinal() ];

    const ptrdiff_t n = j - i ;

    if ( pd.m_size && 0 < n ) {
      array_type::assign(
        (ScalarType*)( (byte_p)(i.m_bucket_ptr->m_entities) +
                       pd.m_base + pd.m_size * i.m_current_entity ),
        (typename array_type::size_type) n );
    }
  }

#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------
/** \brief  \ref stk::mesh::Field "Field" data \ref shards::Array "Array"
 *          for a given array field and bucket
 */
template< typename ScalarType ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct BucketArray< Field<ScalarType,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
  : public shards::ArrayAppend<
  shards::Array<ScalarType,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> ,
  EntityDimension >::type
{
private:
#ifndef DOXYGEN_COMPILE
  typedef unsigned char * byte_p ;
  BucketArray();
  BucketArray( const BucketArray & );
  BucketArray & operator = ( const BucketArray & );
public:

  typedef Field<ScalarType,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> field_type ;

  typedef typename shards::ArrayAppend<
    shards::Array<ScalarType,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> ,
    EntityDimension >::type array_type ;

  BucketArray( const field_type & f , const Bucket & b )
  {
    const Bucket::DataMap & pd = b.m_field_map[ f.mesh_meta_data_ordinal() ];

    if ( pd.m_size ) {
      array_type::assign_stride(
        (ScalarType*)( ((byte_p) b.m_entities ) + pd.m_base ),
        pd.m_stride , (typename array_type::size_type) b.size() );
    }
  }

  BucketArray( const field_type & f,
               const Bucket::iterator & i,
               const Bucket::iterator & j) 
  {
    const ptrdiff_t distance = j - i ;

    if ( 0 < distance ) {

      const Bucket          & b  = * i.m_bucket_ptr ;
      const Bucket::DataMap & pd = b.m_field_map[ f.mesh_meta_data_ordinal() ];

      if ( pd.m_size ) {
        array_type::assign_stride(
          (ScalarType*)( ((byte_p) b.m_entities ) +
                         pd.m_base + pd.m_size * i.m_current_entity ),
          pd.m_stride , (typename array_type::size_type) distance );
      }
    }
  }


#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_FieldData_hpp */

