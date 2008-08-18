/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <stdlib.h>
#include <stdexcept>
#include <sstream>


#define ARRAY_BOUNDS_CHECKING

#include <Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

void array_stride_tag_natural_tag(
  const int rank , const ArrayDimTag ** const dst ,
                   const ArrayDimTag * const * const src )
{
  for ( int i = 0 ; i < rank ; ++i ) { dst[i] = src[ ( rank - 1 ) - i ]; }
}

void array_stride_tag_fortran_tag(
  const int rank , const ArrayDimTag ** const dst ,
                   const ArrayDimTag * const * const src )
{
  for ( int i = 0 ; i < rank ; ++i ) { dst[i] = src[i] ; }
}

size_t array_stride_size( const int rank , const size_t * const stride )
{
  return 0 < rank ? stride[ rank - 1 ] : 0 ;
}

void array_stride_from_natural_sizes(
  const int rank , size_t * const stride , const size_t * const size )
{
  size_t n = 1 ;
  for ( int i = 0 ; i < rank ; ++i ) { stride[i] = n *= size[(rank-1)-i]; }
}

void array_stride_from_fortran_sizes(
  const int rank , size_t * const stride , const size_t * const size )
{
  size_t n = 1 ;
  for ( int i = 0 ; i < rank ; ++i ) { stride[i] = n *= size[i] ; }
}

void array_stride_to_natural_sizes(
  const int rank , const size_t * const stride , size_t * const size )
{
  if ( 0 < rank ) {
    size[ rank - 1 ] = stride[0] ;
    for ( int i = 1 ; i < rank ; ++i ) {
      size[ ( rank - 1 ) - i ] = stride[i] / stride[i-1] ;
    }
  }
}

void array_stride_to_fortran_sizes(
  const int rank , const size_t * const stride , size_t * const size )
{
  if ( 0 < rank ) {
    size[0] = stride[0] ;
    for ( int i = 1 ; i < rank ; ++i ) {
      size[i] = stride[i] / stride[i-1] ;
    }
  }
}

void array_stride_to_natural_index(
  const int rank , const size_t * const stride ,
  const int offset , int * const index )
{
  if ( 0 < rank ) {
    int tmp = offset ;
    for ( int i = rank - 1 ; 0 < i ; ) {
      index[ ( rank - 1 ) - i ] = tmp / stride[i-1] ;
      tmp %= stride[i-1] ;
    }
    index[ rank - 1 ] = tmp ;
  }
}

void array_stride_to_fortran_index(
  const int rank , const size_t * const stride ,
  const int offset , int * const index )
{
  if ( 0 < rank ) {
    int tmp = offset ;
    for ( int i = rank - 1 ; 0 < i ; ) {
      index[i] = tmp / stride[i-1] ;
      tmp %= stride[i-1] ;
    }
    index[0] = tmp ;
  }
}

//----------------------------------------------------------------------

ArrayDimTag::~ArrayDimTag() {}

std::string ArrayDimTag::to_string( size_t n , int i ) const
{
  array_bounds_checking( n , i );
  std::ostringstream s ;
  s << i ;
  return s.str();
}

int ArrayDimTag::to_index( size_t n , const std::string & s ) const
{
  const int i = atoi( s.c_str() );
  array_bounds_checking( n , i );
  return i ;
}

void array_bounds_checking( int arg_size , int arg_ord )
{
  if ( arg_ord < 0 || arg_size <= arg_ord ) {
    std::ostringstream msg ;
    msg << "ARRAY BOUNDS ERROR: Array[" << arg_size ;
    msg << "] ordinal error " << arg_ord ;
    msg << " ; TO DEBUG SET BREAKPOINT IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

void array_bounds_checking( const bool arg_natural ,
                            const int arg_rank ,
                            const size_t * const arg_stride ,
                            const int arg_i1 ,
                            const int arg_i2 ,
                            const int arg_i3 ,
                            const int arg_i4 ,
                            const int arg_i5 ,
                            const int arg_i6 ,
                            const int arg_i7 ,
                            const int arg_i8 )
{
  const int indices[8] =
    { arg_i1 , arg_i2 , arg_i3 , arg_i4 ,
      arg_i5 , arg_i6 , arg_i7 , arg_i8 };

  int sizes[8] ;

  if ( arg_natural ) {
    for ( int ord = 0 ; ord < arg_rank ; ++ord ) {
      const int i = ( arg_rank - 1 ) - ord ;
      sizes[ord] = i ? arg_stride[i] / arg_stride[i-1] : arg_stride[i];
    }
  }
  else {
    for ( int ord = 0 ; ord < arg_rank ; ++ord ) {
      sizes[ord] = ord ? arg_stride[ord] / arg_stride[ord-1] : arg_stride[ord];
    }
  }

  bool ok = true ;

  for ( int ord = 0 ; ok && ord < arg_rank ; ++ord ) {
    ok = 0 <= indices[ord] && indices[ord] < sizes[ord] ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "ARRAY BOUNDS ERROR : Array(" ;
    for ( int ord = 0 ; ord < arg_rank ; ++ord ) {
      if ( ord ) { msg << "," ; }
      msg << sizes[ord] ;
    }
    msg << ") index error (" ;
    for ( int ord = 0 ; ord < arg_rank ; ++ord ) {
      if ( ord ) { msg << "," ; }
      msg << indices[ord] ;
    }
    msg << ")" ;
    msg << " ; TO DEBUG SET BREAKPOINT IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

}

