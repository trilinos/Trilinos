// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdlib.h>
#include <stdexcept>
#include <sstream>

#define SHARDS_ARRAY_BOUNDS_CHECKING

#include <Shards_Array.hpp>

namespace shards {

//----------------------------------------------------------------------

ArrayDimTag::~ArrayDimTag() {}

std::string ArrayDimTag::to_string( ArrayDimTag::size_type n ,
                                    ArrayDimTag::size_type i ) const
{
  array_traits::check_range( i , n );
  std::ostringstream s ;
  s << i ;
  return s.str();
}

ArrayDimTag::size_type
ArrayDimTag::to_index( ArrayDimTag::size_type n , const std::string & s ) const
{
  const int i = atoi( s.c_str() );
  array_traits::check_range( i , n );
  return i ;
}

ArrayDimension::ArrayDimension() {}

ArrayDimension::~ArrayDimension() {}

const char * ArrayDimension::name() const
{ static const char n[] = "anonymous" ; return n ; }

const ArrayDimension & ArrayDimension::tag()
{ static const ArrayDimension self ; return self ; }

//----------------------------------------------------------------------

namespace array_traits {

//----------------------------------------------------------------------

void array_stride_from_fortran_dimensions(
  const array_traits::int_t rank ,
        array_traits::int_t * const stride ,
  const array_traits::int_t * const dim )
{
  array_traits::int_t n = 1 ;
  for ( array_traits::int_t i = 0 ; i < rank ; ++i )
    { stride[i] = n *= dim[i] ; }
}

void array_stride_to_fortran_dimensions(
  const array_traits::int_t rank ,
  const array_traits::int_t * const stride ,
        array_traits::int_t * const dim )
{
  if ( 0 < rank ) {
    dim[0] = stride[0] ;
    for ( array_traits::int_t i = 1 ; i < rank ; ++i ) {
      dim[i] = stride[i] / stride[i-1] ;
    }
  }
}

void array_stride_to_fortran_indices(
  const array_traits::int_t rank ,
  const array_traits::int_t * const stride ,
  const array_traits::int_t offset ,
        array_traits::int_t * const indices )
{
  if ( 0 < rank ) {
    array_traits::int_t tmp = offset ;
    for ( array_traits::int_t i = rank - 1 ; 0 < i ; i--) {
      indices[i] = tmp / stride[i-1] ;
      tmp %= stride[i-1] ;
    }
    indices[0] = tmp ;
  }
}

//----------------------------------------------------------------------

void init_dim(
         int_t dst_stride[] ,
   const int_t src_dimension[] ,
   const int_t rank , const bool natural )
{
  enum { MAX = 8 };
  int_t n = 1 ;
  if ( natural ) {
    for ( int_t k = 0 ; k < rank ; ++k ) {
      dst_stride[k] = ( n *= src_dimension[(rank-1)-k] );
    }
  }
  else {
    for ( int_t k = 0 ; k < rank ; ++k ) {
      dst_stride[k] = ( n *= src_dimension[k] );
    }
  }
  for ( int_t k = rank ; k < MAX ; ++k ) {
    dst_stride[k] = 0 ;
  }
}

void init_tags(
   const ArrayDimTag *       dst_tag[] ,
   const ArrayDimTag * const src_tag[] ,
   const int_t rank , const bool natural )
{
  enum { MAX = 8 };
  if ( natural ) {
    for ( int_t k = 0 ; k < rank ; ++k ) {
      dst_tag[k] = src_tag[(rank-1)-k] ;
    }
  }
  else {
    for ( int_t k = 0 ; k < rank ; ++k ) {
      dst_tag[k] = src_tag[k] ;
    }
  }
  for ( int_t k = rank ; k < MAX ; ++k ) {
    dst_tag[k] = NULL ;
  }
}

void check_rank( const int_t rank ,
                 const int_t test_rank )
{
  if ( rank != test_rank ) {
    std::ostringstream msg ;
    msg << "ARRAY RANK ERROR: ( Rank = " << rank ;
    msg << " ) <= ( test_rank = " << test_rank ;
    msg << " ); THROWN FROM FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

void check_range( const int_t index , const int_t bound )
{
  if ( ! ( 0 <= index && index < bound ) ) {
    std::ostringstream msg ;
    msg << "ARRAY RANGE CHECK FAILED: 0 <= " ;
    msg << index << " < " << bound ;
    msg << " ; THROWN FROM FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

void throw_bad_conversion( const int_t lhs_rank ,
                           const ArrayDimTag * const lhs_tags[] ,
                           const int_t rhs_rank ,
                           const ArrayDimTag * const rhs_tags[] )
{
  std::ostringstream msg ;
  msg << "ARRAY CONVERSION ERROR : Array<T," ;
  msg << lhs_rank ;
  for ( int_t k = 0 ; k < lhs_rank ; ++k ) {
    msg << "," << lhs_tags[k]->name();
  } 
  msg << "> != Array<T," ;
  msg << rhs_rank ;
  for ( int_t k = 0 ; k < lhs_rank ; ++k ) {
    msg << "," << rhs_tags[k]->name();
  } 
  msg << ">" ;
  throw std::runtime_error( msg.str() );
}


void check_indices( const bool arg_natural ,
                    const array_traits::int_t arg_rank ,
                    const array_traits::int_t * const arg_stride ,
                    const array_traits::int_t arg_i1 ,
                    const array_traits::int_t arg_i2 ,
                    const array_traits::int_t arg_i3 ,
                    const array_traits::int_t arg_i4 ,
                    const array_traits::int_t arg_i5 ,
                    const array_traits::int_t arg_i6 ,
                    const array_traits::int_t arg_i7 ,
                    const array_traits::int_t arg_i8 )
{
  const array_traits::int_t indices[8] =
    { arg_i1 , arg_i2 , arg_i3 , arg_i4 ,
      arg_i5 , arg_i6 , arg_i7 , arg_i8 };

  array_traits::int_t sizes[8] ;

  if ( arg_natural ) {
    for ( array_traits::int_t ord = 0 ; ord < arg_rank ; ++ord ) {
      const array_traits::int_t i = ( arg_rank - 1 ) - ord ;
      sizes[ord] = i ? arg_stride[i] / arg_stride[i-1] : arg_stride[i];
    }
  }
  else {
    for ( array_traits::int_t ord = 0 ; ord < arg_rank ; ++ord ) {
      sizes[ord] = ord ? arg_stride[ord] / arg_stride[ord-1] : arg_stride[ord];
    }
  }

  bool ok = true ;

  for ( array_traits::int_t ord = 0 ; ok && ord < arg_rank ; ++ord ) {
    ok = indices[ord] < sizes[ord] ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "ARRAY INDICES ERROR : Array(" ;

    for ( array_traits::int_t ord = 0 ; ord < arg_rank ; ++ord ) {
      if ( ord ) { msg << "," ; }
      if ( sizes[ord] <= indices[ord] ) {
        msg << "dim(" << sizes[ord] << ")<=" ;
      }
      msg << indices[ord] ;
    }
    msg << ")" ;
    msg << " ; THROWN FROM FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

} // namespace array_traits
} // namespace shards

