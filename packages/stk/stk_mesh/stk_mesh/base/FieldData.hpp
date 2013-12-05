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
#include <stk_mesh/base/BulkData.hpp>

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

/** \brief  Size, in bytes, of the field data for each entity */
inline
unsigned field_data_size_per_entity( const FieldBase & f , const Bucket & k )
{
  return BulkData::get(k).field_data_size_per_entity(f, k);
}

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
template< class field_type >
inline
typename FieldTraits< field_type >::data_type *
field_data( const field_type & f , const Entity e )
{
  return stk::mesh::BulkData::get(e).field_data(f, e);
}
#endif

//----------------------------------------------------------------------


#ifndef DOXYGEN_COMPILE

template< typename ScalarType >
struct EntityArray< Field<ScalarType,void,void,void,void,void,void,void> >
  : public shards::Array<ScalarType,shards::RankZero,void,void,void,void,void,void,void>
{
  typedef Field<ScalarType,void,void,void,void,void,void,void> field_type ;
  typedef shards::Array<ScalarType,shards::RankZero,void,void,void,void,void,void,void>
  array_type ;

  EntityArray( const field_type & f , const Bucket& b, Bucket::size_type bucket_ordinal )
    : array_type( BulkData::get(b).field_data( f , b, bucket_ordinal ) ) {}

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

  EntityArray( const field_type & f , const Bucket& b, Bucket::size_type bucket_ordinal) : array_type()
  {
    if ( BulkData::get(b).field_data_size_per_entity(f, b)) {
      array_type::assign_stride( BulkData::get(b).field_data(f, b, bucket_ordinal),
                                 BulkData::get(b).field_data_stride(f, b) );
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
    if (BulkData::get(k).field_data_size_per_entity(f, k)) {
      array_type::assign( BulkData::get(k).field_data(f, k, 0) ,
                          k.size() );
    }
  }

  BucketArray( const field_type & f,
               const Bucket& b,
               const Bucket::iterator i,
               const Bucket::iterator j)
  {
    const ptrdiff_t n = j - i ;

    if ( BulkData::get(b).field_data_size_per_entity(f, b) && 0 < n ) {
      unsigned b_ord = i-b.begin();
      array_type::assign( BulkData::get(b).field_data( f, b, b_ord ),
                          static_cast<typename array_type::size_type>(n) );
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
    if ( BulkData::get(b).field_data_size_per_entity(f, b) ) {
      array_type::assign_stride( BulkData::get(b).field_data(f, b, 0),
                                 BulkData::get(b).field_data_stride(f, b),
                                 static_cast<typename array_type::size_type>(b.size()) );
    }
  }

  BucketArray( const field_type & f,
               const Bucket& b,
               const Bucket::iterator i,
               const Bucket::iterator j)
  {
    const ptrdiff_t distance = j - i ;

    if ( 0 < distance ) {

      if ( BulkData::get(b).field_data_size_per_entity(f, b) ) {
        Bucket::size_type bucket_ordinal = i - b.begin();
        array_type::assign_stride( BulkData::get(b).field_data(f, b, bucket_ordinal),
                                   BulkData::get(b).field_data_stride(f, b),
                                   static_cast<typename array_type::size_type>(distance) );
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

