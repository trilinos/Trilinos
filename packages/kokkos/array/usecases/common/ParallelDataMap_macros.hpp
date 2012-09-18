/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

namespace KokkosArray {

//----------------------------------------------------------------------------

template< typename ValueType , unsigned N1 >
struct PackArray< View< ValueType[][N1] , KOKKOSARRAY_MACRO_DEVICE > , void >
{
  typedef KOKKOSARRAY_MACRO_DEVICE                         device_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE::size_type              size_type ;
  typedef View< ValueType[][N1] , device_type >       array_type ;
  typedef View< ValueType[] , device_type >           buffer_type ;

private:

  buffer_type  output ;
  array_type   input ;
  size_type    base ;

public:

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i ) const
  {
    for ( size_type j = 0 , k = i * N1 ; j < N1 ; ++j , ++k ) {
      output[k] = input(base+i,j);
    }
  }

  inline static
  void pack( const buffer_type & arg_output ,
             const size_type     arg_begin ,
             const size_type     arg_count ,
             const array_type  & arg_input )
  {
    if ( arg_count ) {
      PackArray op ;
      op.output = arg_output ;
      op.input  = arg_input ;
      op.base   = arg_begin ;
      parallel_for( arg_count , op );
    }
  }
};

template< typename ValueType , unsigned N1 >
struct UnpackArray< View< ValueType[][N1] , KOKKOSARRAY_MACRO_DEVICE > , void >
{
  typedef KOKKOSARRAY_MACRO_DEVICE                         device_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE::size_type              size_type ;
  typedef View< ValueType[] , device_type >           buffer_type ;
  typedef View< ValueType[][N1] , device_type >       array_type ;

private:

  array_type   output ;
  buffer_type  input ;
  size_type    base ;

public:

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i ) const
  {
    for ( size_type j = 0 , k = i * N1 ; j < N1 ; ++j , ++k ) {
      output(base+i,j) = input(k);
    }
  }

  inline
  static
  void unpack( const array_type  & arg_output ,
               const buffer_type & arg_input ,
               const size_type     arg_begin ,
               const size_type     arg_count )
  {
    if ( arg_count ) {
      UnpackArray op ;
      op.output = arg_output ;
      op.input  = arg_input ;
      op.base   = arg_begin ;
      parallel_for( arg_count , op );
    }
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct PackArray< View< ValueType[] , KOKKOSARRAY_MACRO_DEVICE > , void >
{
  typedef KOKKOSARRAY_MACRO_DEVICE                device_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE::size_type     size_type ;
  typedef View< ValueType[] , device_type >  array_type ;
  typedef View< ValueType[] , device_type >  buffer_type ;

private:

  buffer_type  output ;
  array_type   input ;
  size_type    base ;

public:

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i ) const
  { output[i] = input(base+i); }

  inline
  static
  void pack( const buffer_type & arg_output ,
             const size_type     arg_begin ,
             const size_type     arg_count ,
             const array_type  & arg_input )
  {
    PackArray op ;
    op.output = arg_output ;
    op.input  = arg_input ;
    op.base   = arg_begin ;
    parallel_for( arg_count , op );
  }
};

template< typename ValueType >
struct UnpackArray< View< ValueType[] , KOKKOSARRAY_MACRO_DEVICE > , void >
{
  typedef KOKKOSARRAY_MACRO_DEVICE                device_type ;
  typedef KOKKOSARRAY_MACRO_DEVICE::size_type     size_type ;
  typedef View< ValueType[] , device_type >  array_type ;
  typedef View< ValueType[] , device_type >  buffer_type ;

private:

  array_type   output ;
  buffer_type  input ;
  size_type    base ;

public:

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i ) const
  { output(base+i) = input[i]; }

  inline
  static
  void unpack( const array_type  & arg_output ,
               const buffer_type & arg_input ,
               const size_type     arg_begin ,
               const size_type     arg_count )
  {
    UnpackArray op ;
    op.output = arg_output ;
    op.input  = arg_input ;
    op.base   = arg_begin ;
    parallel_for( arg_count , op );
  }
};

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

