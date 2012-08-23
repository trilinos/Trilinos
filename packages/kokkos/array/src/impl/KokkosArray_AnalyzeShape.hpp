/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_ANALYZESHAPE_HPP
#define KOKKOSARRAY_ANALYZESHAPE_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< unsigned RankDynamic ,
          unsigned Rank ,
          unsigned ValueSize ,
          unsigned ArgN0 ,
          unsigned ArgN1 ,
          unsigned ArgN2 ,
          unsigned ArgN3 ,
          unsigned ArgN4 ,
          unsigned ArgN5 ,
          unsigned ArgN6 ,
          unsigned ArgN7 >
struct ShapeStaticInfo {
  static const unsigned rank_dynamic = RankDynamic ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = ValueSize ;

  static const unsigned N0 = ArgN0 ;
  static const unsigned N1 = ArgN1 ;
  static const unsigned N2 = ArgN2 ;
  static const unsigned N3 = ArgN3 ;
  static const unsigned N4 = ArgN4 ;
  static const unsigned N5 = ArgN5 ;
  static const unsigned N6 = ArgN6 ;
  static const unsigned N7 = ArgN7 ;
};

// Maximum rank == 8

template< unsigned RankDynamic ,
          unsigned ValueSize ,
          unsigned ArgN0 ,
          unsigned ArgN1 ,
          unsigned ArgN2 ,
          unsigned ArgN3 ,
          unsigned ArgN4 ,
          unsigned ArgN5 ,
          unsigned ArgN6 ,
          unsigned ArgN7 >
struct ShapeStaticInfo<RankDynamic,9,ValueSize,
                       ArgN0,ArgN1,ArgN2,ArgN3,
                       ArgN4,ArgN5,ArgN6,ArgN7> {};

//----------------------------------------------------------------------------

template< class T >
struct AnalyzeShape ;

template< class T >
struct AnalyzeShape {
  typedef       T value_type ;
  typedef       T array_type ;
  typedef       T type ;
  typedef const T const_value_type ;
  typedef const T const_array_type ;
  typedef const T const_type ;
  typedef       T non_const_value_type ;
  typedef       T non_const_array_type ;
  typedef       T non_const_type ;

  typedef ShapeStaticInfo< 0, 0, sizeof(T), 0, 0, 0, 0, 0, 0, 0, 0 >  shape ;
};

template< class T >
struct AnalyzeShape< const T > {
private:
  typedef AnalyzeShape<T> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::const_value_type value_type ;
  typedef typename nested::const_array_type array_type ;
  typedef typename nested::const_type       type ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type const_array_type ;
  typedef typename nested::const_type       const_type ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type non_const_array_type ;
  typedef typename nested::non_const_type       non_const_type ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic ,
                           nested_shape::rank ,
                           nested_shape::value_size ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 ,
                           nested_shape::N7 > shape ;
};

template< class T >
struct AnalyzeShape< T[] >
{
private:
  typedef AnalyzeShape<T> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type array_type [] ;
  typedef typename nested::type       type [] ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type const_array_type [] ;
  typedef typename nested::const_type       const_type [] ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type non_const_array_type [] ;
  typedef typename nested::non_const_type       non_const_type [] ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic + 1 ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           0 ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

template< class T >
struct AnalyzeShape< const T[] >
{
private:
  typedef AnalyzeShape< const T > nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type array_type ;
  typedef typename nested::type       type ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type const_array_type [] ;
  typedef typename nested::const_type       const_type [] ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type non_const_array_type [] ;
  typedef typename nested::non_const_type       non_const_type [] ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic + 1 ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           0 ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

template< class T >
struct AnalyzeShape< T * >
{
private:
  typedef AnalyzeShape<T> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type * array_type ;
  typedef typename nested::type       * type ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type * const_array_type ;
  typedef typename nested::const_type       * const_type ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type * non_const_array_type ;
  typedef typename nested::non_const_type       * non_const_type ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic + 1 ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           0 ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

template< class T >
struct AnalyzeShape< const T * >
{
private:
  typedef AnalyzeShape<const T> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type * array_type ;
  typedef typename nested::type       * type ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type * const_array_type ;
  typedef typename nested::const_type       * const_type ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type * non_const_array_type ;
  typedef typename nested::non_const_type       * non_const_type ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic + 1 ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           0 ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

template< class T , unsigned N >
struct AnalyzeShape< T[N] >
{
private:
  typedef AnalyzeShape<T> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type array_type [N] ;
  typedef typename nested::type  type [N] ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type const_array_type [N] ;
  typedef typename nested::const_type       const_type [N] ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type non_const_array_type [N] ;
  typedef typename nested::non_const_type       non_const_type [N] ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           N ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

template< class T , unsigned N >
struct AnalyzeShape< const T[N] >
{
private:
  typedef AnalyzeShape< const T > nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::value_type value_type ;
  typedef typename nested::array_type array_type [N] ;
  typedef typename nested::type  type [N] ;

  typedef typename nested::const_value_type const_value_type ;
  typedef typename nested::const_array_type const_array_type [N] ;
  typedef typename nested::const_type       const_type [N] ;

  typedef typename nested::non_const_value_type non_const_value_type ;
  typedef typename nested::non_const_array_type non_const_array_type [N] ;
  typedef typename nested::non_const_type       non_const_type [N] ;

  typedef ShapeStaticInfo< nested_shape::rank_dynamic ,
                           nested_shape::rank + 1 ,
                           nested_shape::value_size ,
                           N ,
                           nested_shape::N0 ,
                           nested_shape::N1 ,
                           nested_shape::N2 ,
                           nested_shape::N3 ,
                           nested_shape::N4 ,
                           nested_shape::N5 ,
                           nested_shape::N6 > shape ;
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_ANALYZESHAPE_HPP */

