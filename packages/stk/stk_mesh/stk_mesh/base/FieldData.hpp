/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


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
 *
 * This file contains the portion of the Field API that is related to
 * retrieving field data. It defines the EntityArray and BucketArray
 * structs, which represent the interface between shards arrays and
 * field data. The API for retrieving a raw pointer to the data is
 * simpler, but the EntityArray/BucketArray system can be nice if the
 * Field-type is complex because they handle the striding for you.
 *
 * Note that getting field-data on a per-entity basis is not efficient
 * for a large number of entities; use bucket field-data instead.
 *
 * Examples:
 *   - Get raw data for an entity:
 *     double* data = stk::mesh::field_data(field, entity);
 *
 *   - Get raw data for a bucket:
 *     double* data = stk::mesh::field_data(field, bucket);
 *
 *   - Get EntityArray data for an entity, assuming field type: Field<double, Cartesian>:
 *     EntityArray<Field<double, Cartesian> > data_array(field, entity);
 *
 *   - Get BucketArray data for a bucket, assuming field type: Field<double, Cartesian>:
 *     BucketArray<Field<double, Cartesian> > data_array(field, bucket);
 *
 *   - Using (Entity|Bucket)Array, assuming we are dealing with coordinates:
 *       int num_coords_per_node = data_array.dimension(0);
 *       int num_nodes_in_bucket = data_array.dimension(1);
 *       for (int n = 0; n < num_nodes_in_bucket; ++n) {
 *         for (int c = 0; c < num_coords_per_node; ++c) {
 *           cout << data_data(c, n) << ", ";
 *         }
 *         cout << endl;
 *       }
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
  return k.field_data_size(f);
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
  return i.m_bucket_ptr->field_data( f, *i );
}


/** \brief  Pointer to the field data array */
template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Entity & e )
{
  return e.bucket().field_data( f, e );
}

//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

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
    if (b.field_data_size(f)) {
      array_type::assign_stride(
          (ScalarType*)(b.field_data_location(f, e ) ), 
          b.field_data_stride(f)
          );
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
    if (k.field_data_size(f)) {
      array_type::assign( (ScalarType*)( k.field_data_location(f,k[0]) ) ,
                          k.size() );

    }
  }

  BucketArray( const field_type & f,
               const Bucket::iterator & i,
               const Bucket::iterator & j) 
  {
    const ptrdiff_t n = j - i ;

    if ( i.m_bucket_ptr->field_data_size(f) && 0 < n ) {
      array_type::assign(
        (ScalarType*)( i.m_bucket_ptr->field_data_location( f, *i ) ),
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
    if ( b.field_data_size(f) ) {
      array_type::assign_stride(
        (ScalarType*)( b.field_data_location(f,b[0]) ),
        b.field_data_stride(f) , (typename array_type::size_type) b.size() );
    }
  }

  BucketArray( const field_type & f,
               const Bucket::iterator & i,
               const Bucket::iterator & j) 
  {
    const ptrdiff_t distance = j - i ;

    if ( 0 < distance ) {

      const Bucket          & b  = * i.m_bucket_ptr ;

      if ( b.field_data_size(f) ) {
        array_type::assign_stride(
          (ScalarType*)( b.field_data_location(f,*i) ),
          b.field_data_stride(f) , (typename array_type::size_type) distance );
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

