/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United Traitss Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_base_FieldTraits_hpp
#define stk_mesh_base_FieldTraits_hpp

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>

#include <Shards_Array.hpp>

namespace stk {
namespace mesh {

/**
 * FieldTraits provide an API for making queries about field types.
 * Examples:
 *   - Get the scalar data type contained by a field
 *     stk::mesh::FieldTraits< field_type >::data_type
 *   - Get the dimensional rank (number of dimensions) in a field type
 *     stk::mesh::FieldTraits< field_type >::Rank
 */

template<>
struct FieldTraits<FieldBase>
{
public:
  typedef shards::array_traits::Helper<void,shards::RankZero,
                                       void,void,void,void,void,void,void,void>
    Helper ;

  typedef void data_type ; ///< \brief  Data type of the field's members
  typedef void   tag1 ;      ///< \brief  Array dimension tag
  typedef void   tag2 ;      ///< \brief  Array dimension tag
  typedef void   tag3 ;      ///< \brief  Array dimension tag
  typedef void   tag4 ;      ///< \brief  Array dimension tag
  typedef void   tag5 ;      ///< \brief  Array dimension tag
  typedef void   tag6 ;      ///< \brief  Array dimension tag
  typedef void   tag7 ;      ///< \brief  Array dimension tag

  /** \brief  Multidimensional array rank */
  enum { Rank = 0 };

  static void assign_tags( const shards::ArrayDimTag ** tags ) {}
};


/** \brief  Scalar type and multi-dimensional array traits of a Field */
template< typename Scalar >
struct FieldTraits< Field<Scalar,void,void,void,void,void,void,void> >
{
public:
  typedef shards::array_traits::Helper<Scalar,shards::RankZero,
                                       void,void,void,void,void,void,void,void>
    Helper ;

  typedef Scalar data_type ; ///< \brief  Data type of the field's members
  typedef void   tag1 ;      ///< \brief  Array dimension tag
  typedef void   tag2 ;      ///< \brief  Array dimension tag
  typedef void   tag3 ;      ///< \brief  Array dimension tag
  typedef void   tag4 ;      ///< \brief  Array dimension tag
  typedef void   tag5 ;      ///< \brief  Array dimension tag
  typedef void   tag6 ;      ///< \brief  Array dimension tag
  typedef void   tag7 ;      ///< \brief  Array dimension tag

  /** \brief  Multidimensional array rank */
  enum { Rank = 0 };

  static void assign_tags( const shards::ArrayDimTag ** tags ) {}
};

/** \brief  Scalar type and multi-dimensional array traits of a Field */
template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct FieldTraits< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{
public:
  typedef shards::array_traits::Helper<Scalar,shards::FortranOrder,
                                       Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
    Helper ;

  typedef Scalar data_type ; ///< \brief  Data type of the field's members
  typedef Tag1   tag1 ;      ///< \brief  Array dimension tag
  typedef Tag2   tag2 ;      ///< \brief  Array dimension tag
  typedef Tag3   tag3 ;      ///< \brief  Array dimension tag
  typedef Tag4   tag4 ;      ///< \brief  Array dimension tag
  typedef Tag5   tag5 ;      ///< \brief  Array dimension tag
  typedef Tag6   tag6 ;      ///< \brief  Array dimension tag
  typedef Tag7   tag7 ;      ///< \brief  Array dimension tag

  /** \brief  Multidimensional array rank */
  enum { Rank = Helper::Rank };

  static void assign_tags( const shards::ArrayDimTag ** tags )
    { Helper::assign_tags( tags ); }
};


} //namespace mesh
} //namespace stk

#endif //stk_mesh_base_FieldTraits_hpp
