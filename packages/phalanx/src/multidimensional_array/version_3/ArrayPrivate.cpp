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
 */

#include <stdlib.h>
#include <stdexcept>
#include <sstream>

#define ARRAY_BOUNDS_CHECKING

#include <Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

ArrayDimTag::~ArrayDimTag() {}

std::string ArrayDimTag::to_string( unsigned n ,
                                    unsigned i ) const
{
  array_check_index( n , i );
  std::ostringstream s ;
  s << i ;
  return s.str();
}

unsigned ArrayDimTag::to_index( unsigned n ,
                                        const std::string & s ) const
{
  const int i = atoi( s.c_str() );
  array_check_index( n , i );
  return i ;
}

ArrayDimension::~ArrayDimension() {}

const char * ArrayDimension::name() const
{ static const char n[] = "anonymous" ; return n ; }

const ArrayDimension & ArrayDimension::tag()
{ static const ArrayDimension self ; return self ; }

//----------------------------------------------------------------------

void array_stride_from_natural_dimensions(
  const unsigned rank ,
        unsigned * const stride ,
  const unsigned * const dim )
{
  unsigned n = 1 ;
  for ( unsigned i = 0 ; i < rank ; ++i )
    { stride[i] = n *= dim[(rank-1)-i]; }
}

unsigned array_stride_size(
  const unsigned  rank ,
  const unsigned * const stride )
{ return 0 < rank ? stride[ rank - 1 ] : 0 ; }

void array_stride_to_natural_dimensions(
  const unsigned   rank ,
  const unsigned  * const stride ,
        unsigned * const dim )
{
  if ( 0 < rank ) {
    dim[ rank - 1 ] = stride[0] ;
    for ( unsigned i = 1 ; i < rank ; ++i ) {
      dim[ ( rank - 1 ) - i ] = stride[i] / stride[i-1] ;
    }
  }
}

void array_stride_to_natural_indices(
  const unsigned   rank ,
  const unsigned  * const stride ,
  const unsigned    offset ,
        unsigned * const indices )
{
  if ( 0 < rank ) {
    unsigned tmp = offset ;
    for ( unsigned i = rank - 1 ; 0 < i ; ) {
      indices[ ( rank - 1 ) - i ] = tmp / stride[i-1] ;
      tmp %= stride[i-1] ;
    }
    indices[ rank - 1 ] = tmp ;
  }
}

//----------------------------------------------------------------------

void array_stride_from_fortran_dimensions(
  const unsigned rank ,
        unsigned * const stride ,
  const unsigned * const dim )
{
  unsigned n = 1 ;
  for ( unsigned i = 0 ; i < rank ; ++i )
    { stride[i] = n *= dim[i] ; }
}

void array_stride_to_fortran_dimensions(
  const unsigned rank ,
  const unsigned * const stride ,
        unsigned * const dim )
{
  if ( 0 < rank ) {
    dim[0] = stride[0] ;
    for ( unsigned i = 1 ; i < rank ; ++i ) {
      dim[i] = stride[i] / stride[i-1] ;
    }
  }
}

void array_stride_to_fortran_indices(
  const unsigned rank ,
  const unsigned * const stride ,
  const unsigned offset ,
        unsigned * const indices )
{
  if ( 0 < rank ) {
    unsigned tmp = offset ;
    for ( unsigned i = rank - 1 ; 0 < i ; ) {
      indices[i] = tmp / stride[i-1] ;
      tmp %= stride[i-1] ;
    }
    indices[0] = tmp ;
  }
}

//----------------------------------------------------------------------

void array_check_rank( const unsigned rank ,
                       const unsigned test_rank )
{
  if ( rank != test_rank ) {
    std::ostringstream msg ;
    msg << "ARRAY RANK ERROR: ( Rank = " << rank ;
    msg << " ) <= ( test_rank = " << test_rank ;
    msg << " ); THROWN IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

void array_check_ordinal( const unsigned rank ,
                          const unsigned test_ordinal )
{
  if ( rank <= test_ordinal ) {
    std::ostringstream msg ;
    msg << "ARRAY ORDINAL ERROR: ( Rank = " << rank ;
    msg << " ) <= ( ordinal = " << test_ordinal ;
    msg << " ); THROWN IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

void array_check_index( const unsigned dim ,
                        const unsigned index )
{
  if ( dim <= index ) {
    std::ostringstream msg ;
    msg << "ARRAY INDEX ERROR: ( Dim = " << dim ;
    msg << " ) <= ( index = " << index ;
    msg << " ); THROWN IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

void array_check_offset( const unsigned size ,
                         const unsigned test_offset )
{
  if ( size <= test_offset ) {
    std::ostringstream msg ;
    msg << "ARRAY OFFSET ERROR: ( Size = " << size ;
    msg << " ) <= ( ordinal = " << test_offset ;
    msg << " ); THROWN IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

void array_check_indices( const bool arg_natural ,
                          const unsigned arg_rank ,
                          const unsigned * const arg_stride ,
                          const unsigned arg_i1 ,
                          const unsigned arg_i2 ,
                          const unsigned arg_i3 ,
                          const unsigned arg_i4 ,
                          const unsigned arg_i5 ,
                          const unsigned arg_i6 ,
                          const unsigned arg_i7 ,
                          const unsigned arg_i8 )
{
  const unsigned indices[8] =
    { arg_i1 , arg_i2 , arg_i3 , arg_i4 ,
      arg_i5 , arg_i6 , arg_i7 , arg_i8 };

  unsigned sizes[8] ;

  if ( arg_natural ) {
    for ( unsigned ord = 0 ; ord < arg_rank ; ++ord ) {
      const unsigned i = ( arg_rank - 1 ) - ord ;
      sizes[ord] = i ? arg_stride[i] / arg_stride[i-1] : arg_stride[i];
    }
  }
  else {
    for ( unsigned ord = 0 ; ord < arg_rank ; ++ord ) {
      sizes[ord] = ord ? arg_stride[ord] / arg_stride[ord-1] : arg_stride[ord];
    }
  }

  bool ok = true ;

  for ( unsigned ord = 0 ; ok && ord < arg_rank ; ++ord ) {
    ok = indices[ord] < sizes[ord] ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "ARRAY INDICES ERROR : Array(" ;
    for ( unsigned ord = 0 ; ord < arg_rank ; ++ord ) {
      if ( ord ) { msg << "," ; }
      msg << sizes[ord] ;
    }
    msg << ") given index (" ;
    for ( unsigned ord = 0 ; ord < arg_rank ; ++ord ) {
      if ( ord ) { msg << "," ; }
      msg << indices[ord] ;
    }
    msg << ")" ;
    msg << " ; THROWN IN FILE " ;
    msg << __FILE__ << " AT LINE " << __LINE__ ;
    throw std::runtime_error( msg.str() );
  }
}

}

