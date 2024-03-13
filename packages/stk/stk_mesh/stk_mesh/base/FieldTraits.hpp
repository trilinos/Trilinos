// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
struct STK_DEPRECATED FieldTraits<FieldBase>
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
struct STK_DEPRECATED FieldTraits< Field<Scalar,void,void,void,void,void,void,void> >
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
struct STK_DEPRECATED FieldTraits< Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
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
