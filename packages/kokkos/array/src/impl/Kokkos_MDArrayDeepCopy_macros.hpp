/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_FUNCTION)         || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including impl/Kokkos_MDArrayDeepCopy_macros.hpp without macros defined"

#else

#include <Kokkos_ParallelFor.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class ValueType , class DeviceType ,
          class MapDst , class MapSrc , unsigned Rank >
class MDArrayDeepCopyMember ;

template< class ValueType , class DeviceType ,
          class MapDst , class MapSrc , unsigned Rank >
class MDArrayDeepCopyFunctor ;

//----------------------------------------------------------------------------

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,8>
{
public:
  enum { RANK = 8 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork ,
                    const dst_type & dst ,
                    const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] , indices[6] , indices[7] ) =

    src( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] , indices[6] , indices[7] );

  }
};

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,7>
{
public:
  enum { RANK = 7 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] , indices[6] ) =

    src( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] , indices[6] );

  }
};

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,6>
{
public:
  enum { RANK = 6 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] ) =

    src( indices[0] , indices[1] , indices[2] , indices[3] ,
         indices[4] , indices[5] );

  }
};

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,5>
{
public:
  enum { RANK = 5 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] , indices[3] , indices[4] ) =
    src( indices[0] , indices[1] , indices[2] , indices[3] , indices[4] );
  }
};


template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,4>
{
public:
  enum { RANK = 4 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] , indices[3] ) =
    src( indices[0] , indices[1] , indices[2] , indices[3] );
  }
};

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,3>
{
public:
  enum { RANK = 3 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] , indices[2] ) =
    src( indices[0] , indices[1] , indices[2] );
  }
};

template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,2>
{
public:
  enum { RANK = 2 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  {
    size_type indices[ RANK ];

    dst.inverse_map( iwork , indices );

    dst( indices[0] , indices[1] ) =
    src( indices[0] , indices[1] );
  }
};


template< class ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopyMember<ValueType,KOKKOS_MACRO_DEVICE,MapDst,MapSrc,1>
{
public:
  enum { RANK = 1 };

  typedef KOKKOS_MACRO_DEVICE :: size_type size_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst_type ;
  typedef MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src_type ;

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void copy( size_type iwork , const dst_type & dst , const src_type & src )
  { dst( iwork ) = src( iwork ); }
};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class ValueType , class MapDst , class MapSrc , unsigned Rank >
class MDArrayDeepCopyFunctor< ValueType , KOKKOS_MACRO_DEVICE , MapDst , MapSrc , Rank >
{
private:
  enum { OK = StaticAssert< ! SameType<MapDst,MapSrc>::value >::value };
public :
  typedef KOKKOS_MACRO_DEVICE                      device_type ;
  typedef typename device_type::size_type size_type ;
  typedef MDArrayDeepCopyMember< ValueType , KOKKOS_MACRO_DEVICE , MapDst , MapSrc , Rank > member ;

  MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > dst ;
  MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > src ;

  MDArrayDeepCopyFunctor( MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapDst > arg_dst ,
                          MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapSrc > arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { member::copy( iwork , dst , src ); }
};

//----------------------------------------------------------------------------

/** \brief  Same map and contiguous memory */
template< class ValueType , class MapOption >
class MDArrayDeepCopyFunctor< ValueType , KOKKOS_MACRO_DEVICE ,
                              MapOption , MapOption , 0 > {
public:

  typedef KOKKOS_MACRO_DEVICE                      device_type ;
  typedef typename device_type::size_type size_type ;

        ValueType * const dst_ptr ;
  const ValueType * const src_ptr ;

  MDArrayDeepCopyFunctor( MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapOption > arg_dst ,
                          MDArrayView< ValueType , KOKKOS_MACRO_DEVICE , MapOption > arg_src )
    : dst_ptr( arg_dst.ptr_on_device() ), src_ptr( arg_src.ptr_on_device() ) {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { dst_ptr[ iwork ] = src_ptr[ iwork ]; }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_HOSTMDARRAYDEEPCOPY_HPP */


