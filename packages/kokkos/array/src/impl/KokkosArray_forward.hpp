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

#ifndef KOKKOS_IMPL_FORWARD_HPP
#define KOKKOS_IMPL_FORWARD_HPP

//----------------------------------------------------------------------------

namespace KokkosArray {

class Host ;

template < class Device > struct HostMapped ;

} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class A , class B >
struct SameType { static const bool value = false ; };

template< class A >
struct SameType<A,A> { static const bool value = true ; };

struct MirrorUseView {};

template< class DstType , class SrcType >  struct Factory ;

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType > class MMultiply ;
template < class ValueType, class Device > class MatrixMarketWriter ;
template < class MultiVectorType > class Update ;

}
}

//----------------------------------------------------------------------------

namespace KokkosArray {

template< class ArrayType >
inline
typename ArrayType::HostMirror
create_mirror( const ArrayType & input , const Impl::MirrorUseView )
{
  typedef typename ArrayType::HostMirror MirrorType ;
  return Impl::Factory< MirrorType , Impl::MirrorUseView >::create( input );
}

template< class ArrayType >
inline
typename ArrayType::HostMirror
create_mirror( const ArrayType & input )
{
  typedef typename ArrayType::HostMirror MirrorType ;
#if KOKKOS_MIRROR_VIEW_OPTIMIZE
  return Impl::Factory< MirrorType , Impl::MirrorUseView >::create( input );
#else
  return Impl::Factory< MirrorType , ArrayType >::create( input );
#endif
}

}

#if 0

namespace KokkosArray {

//----------------------------------------------------------------------------
// Possible consolidation of 'create_*' function APIs

template< class ArrayType >
inline
ArrayType create()
{ return Impl::Factory<ArrayType,void>::create(); }

template< class ArrayType , class Arg0Type>
inline
ArrayType create( const Arg0Type & arg0 )
{ return Impl::Factory<ArrayType,void>::create( arg0 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 )
{ return Impl::Factory<ArrayType,void>::create( arg0 , arg1 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 )
{ return Impl::Factory<ArrayType,void>::create( arg0 , arg1 , arg2 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 )
{ return Impl::Factory<ArrayType,void>::create( arg0 , arg1 , arg2 , arg3 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type , class Arg4Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 ,
                  const Arg4Type & arg4 )
{ return Impl::Factory<ArrayType,void>
             ::create( arg0 , arg1 , arg2 , arg3 , arg4 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type , class Arg4Type , class Arg5Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 ,
                  const Arg4Type & arg4 , const Arg5Type & arg5 )
{ return Impl::Factory<ArrayType,void>
             ::create( arg0 , arg1 , arg2 , arg3 , arg4 , arg5 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type , class Arg4Type , class Arg5Type ,
          class Arg6Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 ,
                  const Arg4Type & arg4 , const Arg5Type & arg5 ,
                  const Arg4Type & arg6 )
{ return Impl::Factory<ArrayType,void>
             ::create( arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 ); }

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type , class Arg4Type , class Arg5Type ,
          class Arg6Type , class Arg7Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 ,
                  const Arg4Type & arg4 , const Arg5Type & arg5 ,
                  const Arg4Type & arg6 , const Arg7Type & arg7 )
{
  return Impl::Factory<ArrayType,void>
             ::create( arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );
}

template< class ArrayType ,
          class Arg0Type , class Arg1Type , class Arg2Type ,
          class Arg3Type , class Arg4Type , class Arg5Type ,
          class Arg6Type , class Arg7Type , class Arg8Type >
inline
ArrayType create( const Arg0Type & arg0 , const Arg1Type & arg1 ,
                  const Arg2Type & arg2 , const Arg3Type & arg3 ,
                  const Arg4Type & arg4 , const Arg5Type & arg5 ,
                  const Arg4Type & arg6 , const Arg7Type & arg7 ,
                  const Arg4Type & arg8 )
{
  return Impl::Factory<ArrayType,void>
             ::create( arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8 );
}

} // namespace KokkosArray

#endif /* possible consolidation of 'create_*' APIs */


//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_FORWARD_HPP */



