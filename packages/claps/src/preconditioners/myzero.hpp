//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef MYZERO_HPP
#define MYZERO_HPP
#include <string.h>

inline void myzero( int *data, size_t num )
{
  memset( data, 0, num*sizeof(int) );
}

inline void myzero(short *data, size_t num )
{
  memset( data, 0, num*sizeof(short) );
}

inline void myzero(unsigned short *data, size_t num )
{
  memset( data, 0, num*sizeof(unsigned short) );
}

inline void myzero( float *data, size_t num )
{
  memset( data, 0, num*sizeof(float) );
}

inline void myzero( double *data, size_t num )
{
  memset( data, 0, num*sizeof(double) );
}

inline void myzero( signed char *data, size_t num )
{
  memset( data, 0, num*sizeof(signed char) );
}

inline void myzero( unsigned char *data, size_t num )
{
  memset( data, 0, num*sizeof(unsigned char) );
}
#endif // MYZERO_HPP
