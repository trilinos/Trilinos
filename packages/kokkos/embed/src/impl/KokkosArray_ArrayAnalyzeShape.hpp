/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_ARRAYANALYZESHAPE_HPP
#define KOKKOSARRAY_ARRAYANALYZESHAPE_HPP

#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_AnalyzeShape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  The type 'T' is expected be an intrinsic scalar type;
 *          however, the nesting analysis is invoked for extensibility
 *          to other embedded types.  If not an intrinsic scalar type
 *          then be cautious of side effects.
 */

template< typename T , unsigned N >
struct AnalyzeShape< const Array< T , N > >
  : public ShapeInsert< typename AnalyzeShape< const T >::shape , N >::type
{
private:
  typedef AnalyzeShape< const T > nested ;
public:

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::scalar_type  scalar_type ;
  typedef typename nested::array_type   array_type[ N ];
  typedef const Array< T , N , void >   value_type ;
  typedef const Array< T , N , void >   type ;

  typedef scalar_type const_scalar_type ;
  typedef array_type  const_array_type ;
  typedef value_type  const_value_type ;
  typedef type        const_type ;

  typedef typename nested::non_const_scalar_type   non_const_scalar_type ;
  typedef typename nested::non_const_array_type    non_const_array_type[ N ];
  typedef Array< T , N , void >                    non_const_value_type ;
  typedef Array< T , N , void >                    non_const_type ;
};

template< typename T , unsigned N >
struct AnalyzeShape< Array< T , N > >
  : public ShapeInsert< typename AnalyzeShape<T>::shape , N >::type
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::scalar_type    scalar_type ;
  typedef typename nested::array_type     array_type[ N ];
  typedef          Array< T , N , void >  value_type ;
  typedef          Array< T , N , void >  type ;

  typedef typename nested::const_scalar_type    const_scalar_type ;
  typedef typename nested::const_array_type     const_array_type[ N ];
  typedef          const Array< T , N , void >  const_value_type ;
  typedef          const Array< T , N , void >  const_type ;

  typedef          scalar_type  non_const_scalar_type ;
  typedef          array_type   non_const_array_type ;
  typedef          value_type   non_const_value_type ;
  typedef          type         non_const_type ;
};

//----------------------------------------------------------------------------

template< typename T >
struct AnalyzeShape< const Array< T , 0 > >
  : public ShapeInsert< typename AnalyzeShape< const T >::shape , 0 >::type
{
private:
  typedef AnalyzeShape< const T > nested ;
public:

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::scalar_type    scalar_type ;
  typedef typename nested::array_type   * array_type ;
  typedef const Array< T , 0 , void >     value_type ;
  typedef const Array< T , 0 , void >     type ;

  typedef scalar_type const_scalar_type ;
  typedef array_type  const_array_type ;
  typedef value_type  const_value_type ;
  typedef type        const_type ;

  typedef typename nested::non_const_scalar_type   non_const_scalar_type ;
  typedef typename nested::non_const_array_type  * non_const_array_type ;
  typedef Array< T , 0 , void >                    non_const_value_type ;
  typedef Array< T , 0 , void >                    non_const_type ;
};

template< typename T >
struct AnalyzeShape< Array< T , 0 > >
  : public ShapeInsert< typename AnalyzeShape<T>::shape , 0 >::type
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::scalar_type    scalar_type ;
  typedef typename nested::array_type   * array_type ;
  typedef          Array< T , 0 , void >  value_type ;
  typedef          Array< T , 0 , void >  type ;

  typedef typename nested::const_scalar_type    const_scalar_type ;
  typedef typename nested::const_array_type   * const_array_type ;
  typedef          const Array< T , 0 , void >  const_value_type ;
  typedef          const Array< T , 0 , void >  const_type ;

  typedef          scalar_type  non_const_scalar_type ;
  typedef          array_type   non_const_array_type ;
  typedef          value_type   non_const_value_type ;
  typedef          type         non_const_type ;
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_ARRAYANALYZESHAPE_HPP */

