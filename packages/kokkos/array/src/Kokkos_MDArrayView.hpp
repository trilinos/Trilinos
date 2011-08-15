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

#ifndef KOKKOS_MDARRAYVIEW_HPP
#define KOKKOS_MDARRAYVIEW_HPP

#include <cstddef>
#include <string>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>
#include <impl/Kokkos_MDArrayIndexMap.hpp>
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>

namespace Kokkos {

enum { MDArrayMaxRank = 8 };

template< typename ValueType , class DeviceType >
class MDArrayView ;

template< typename ValueType , class DeviceType >
MDArrayView< ValueType , DeviceType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 );

namespace Impl {

template< typename ValueType , class DeviceDst , class DeviceSrc ,
          bool same_memory_space =
              SameType< typename DeviceDst::memory_space ,
                        typename DeviceSrc::memory_space >::value ,
          bool same_mdarray_map =
              SameType< typename DeviceDst::mdarray_map ,
                        typename DeviceSrc::mdarray_map >::value ,
          bool both_contiguous =
              MDArrayView< ValueType , DeviceDst >::Contiguous &&
              MDArrayView< ValueType , DeviceSrc >::Contiguous >
class MDArrayDeepCopy ;

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Multidimensional array allocated and mapped
 *          onto the memory space of a compute device.
 *
 *  The array is a simple rank-N container of simple scalar values
 *  where 1 <= N <= 8.
 *
 *  The first rank is the parallel work index.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped array and thus achieve portability
 *  across compute devices.
 *
 *  Several APIs for MDArray creation functions are available.
 *  The "labeled" group creates MDArray's with string labels that
 *  will appear in error messages.
 *
 *  create_labeled_mdarray< ValueType , DeviceType >( label , nP , ... );
 *  create_labeled_mdarray< MDArrayView<...> >( label , nP , ... );
 *
 *  The "unlabeled" group creates MDArray's with NULL string labels.
 *
 *  create_mdarray< ValueType , DeviceType >( nP , ... );
 *  create_mdarray< MDArrayView<...> >( nP , ... );
 */

template< typename ValueType , class DeviceType = DeviceHost >
class MDArrayView {
public:
  typedef ValueType                          value_type ;
  typedef DeviceType                         device_type ;
  typedef typename DeviceType::mdarray_map   mdarray_map ;
  typedef typename DeviceType::memory_space  memory_space ;
  typedef typename DeviceType::size_type     size_type ;

  typedef MDArrayView< value_type ,
                       Serial< HostMemory , mdarray_map > > HostView ;

  /*------------------------------------------------------------------*/
  /** \brief  True if the array type has contiguous memory
   *          If contiguous then can get a pointer to the memory.
   */
  enum { Contiguous =
         Impl::MDArrayIndexMap< memory_space , mdarray_map >::Contiguous };

  /** \brief  Query rank of the array */
  size_type rank() const { return m_map.rank(); }

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  size_type dimension( const iType & rank_ordinate ) const
    { return m_map.dimension( rank_ordinate ); }

  /** \brief  Query all dimensions */
  template< typename iType >
  void dimensions( iType * const dims ) const
    { return m_map.dimensions( dims ); }

  /** \brief  Query total number of members */
  size_type size() const { return m_map.size(); }

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5,i6,i7) ]; }

  /** \brief  Query value of a rank 7 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5,i6) ]; }

  /** \brief  Query value of a rank 6 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const 
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4,i5) ]; }

  /** \brief  Query value of a rank 5 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3,i4) ]; }

  /** \brief  Query value of a rank 4 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2,i3) ]; }

  /** \brief  Query value of a rank 3 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1,i2) ]; }

  /** \brief  Query value of a rank 2 array */
  template< typename iTypeP , typename iType1 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP,i1) ]; }

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  value_type & operator()( const iTypeP & iP ) const
  { return m_memory.ptr_on_device()[ m_map.offset(iP) ]; }

  /*------------------------------------------------------------------*/
  /** \brief  Memory is contiguous, OK to return pointer */
  inline
  value_type * ptr_on_device() const { return m_memory.ptr_on_device(); }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArrayView() : m_memory(), m_map() {}

  /** \brief  Construct another view of the 'rhs' array */
  inline
  MDArrayView( const MDArrayView & rhs )
    : m_memory() , m_map( rhs.m_map )
    { memory_space::assign_memory_view( m_memory , rhs.m_memory ); }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  MDArrayView & operator = ( const MDArrayView & rhs )
    {
      memory_space::assign_memory_view( m_memory , rhs.m_memory );
      m_map = rhs.m_map ;
      return *this ;
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  ~MDArrayView() {}

  /*------------------------------------------------------------------*/
  /** \brief  Query if non-NULL view */
  operator bool () const 
  { return m_memory.operator bool(); }

  bool operator == ( const MDArrayView & rhs ) const
  { return m_memory.operator == ( rhs.m_memory ); }

  bool operator != ( const MDArrayView & rhs ) const
  { return m_memory.operator != ( rhs.m_memory ); }

private:
  enum { OK_memory_space =
           Impl::StaticAssert<
             Impl::SameType< HostMemory , memory_space >::value >::value };


  MemoryView< value_type , memory_space >            m_memory ;
  Impl::MDArrayIndexMap< memory_space , mdarray_map > m_map ;

  MDArrayView( const std::string & label ,
               size_t nP , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 )
    : m_memory()
    , m_map( nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 )
    {
      memory_space::allocate_memory_view( m_memory , m_map.size() , label );
      value_type * dst = m_memory.ptr_on_device();
      value_type * const dst_end = dst + m_map.size();
      while ( dst_end != dst ) *dst++ = 0 ;
    }

  template< typename V , class D >
  friend
  MDArrayView< V , D >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class DeviceDst , class DeviceSrc ,
            bool , bool , bool >
  friend
  class Impl::MDArrayDeepCopy ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  This is THE creation function.
 *          All other versions call this function.
 */
template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                        size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return MDArrayView< ValueType , DeviceType >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------
/** \brief  MDArrayView compatibility traits. */

template< class MDArrayViewLHS , class MDArrayViewRHS >
class MDArrayViewTraits {
public:

  /** \brief  Does their data reside in the same memory space ? */
  enum { same_memory_space =
           Impl::SameType< typename MDArrayViewLHS::memory_space ,
                           typename MDArrayViewRHS::memory_space >::value };

  /** \brief  Does their data have the same mapping ? */
  enum { same_mdarray_map =
           Impl::SameType< typename MDArrayViewLHS::mdarray_map ,
                           typename MDArrayViewRHS::mdarray_map >::value };

  /** \brief  Do they have the same value type ? */
  enum { same_value_type =
           Impl::SameType< typename MDArrayViewLHS::value_type ,
                           typename MDArrayViewRHS::value_type >::value };

  /** \brief  Are the types compatible for an assignment operation? */
  enum { assignment_available = same_value_type };

  /** \brief  Does assignment have view semantics (shallow copy)? */
  enum { assignment_view = same_memory_space && same_mdarray_map };

  /** \brief  Does assignment require a deep_copy? */
  enum { assignment_copy = ! same_memory_space || ! same_mdarray_map };

  /** \brief  Does assignment require a permutation of data? */
  enum { assignment_copy_permutation = ! same_mdarray_map };

  /** \brief  Are the dimensions equal? */
  inline
  static bool equal_dimension( const MDArrayViewLHS & lhs ,
                               const MDArrayViewRHS & rhs )
    {
      typedef typename MDArrayViewLHS::size_type size_type ;

      return lhs.rank()       == (size_type) rhs.rank() &&
             lhs.dimension(0) == (size_type) rhs.dimension(0) &&
             lhs.dimension(1) == (size_type) rhs.dimension(1) &&
             lhs.dimension(2) == (size_type) rhs.dimension(2) &&
             lhs.dimension(3) == (size_type) rhs.dimension(3) &&
             lhs.dimension(4) == (size_type) rhs.dimension(4) &&
             lhs.dimension(5) == (size_type) rhs.dimension(5) &&
             lhs.dimension(6) == (size_type) rhs.dimension(6) &&
             lhs.dimension(7) == (size_type) rhs.dimension(7);
    }
};

template< class MDArrayViewLHS , class MDArrayViewRHS >
inline
MDArrayViewTraits< MDArrayViewLHS , MDArrayViewRHS >
mdarray_view_traits( const MDArrayViewLHS & , const MDArrayViewRHS & )
{ return MDArrayViewTraits< MDArrayViewLHS , MDArrayViewRHS >(); }

template< class MDArrayViewLHS , class MDArrayViewRHS >
inline
MDArrayViewTraits< MDArrayViewLHS , MDArrayViewRHS >
mdarray_view_equal_dimension( const MDArrayViewLHS & lhs , const MDArrayViewRHS & rhs )
{ return MDArrayViewTraits< MDArrayViewLHS , MDArrayViewRHS >::equal_dimension(lhs,rhs); }

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class DeviceSrc >
inline
void deep_copy( const MDArrayView<ValueType,DeviceDst> & dst ,
                const MDArrayView<ValueType,DeviceSrc> & src )
{
  typedef MDArrayView<ValueType,DeviceDst>  dst_type ;
  typedef MDArrayView<ValueType,DeviceSrc>  src_type ;

  Impl::mdarray_require_equal_dimension( dst , src );

  Impl::MDArrayDeepCopy< ValueType, DeviceDst, DeviceSrc >::run( dst , src );
}

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
  
template< class DeviceType , class MDArrayDst ,
                             class MDArraySrc , unsigned Rank >
class MDArrayDeepCopyFunctor ;

//----------------------------------------------------------------------------

/** \brief  Fully conformal and contiguous arrays */
template< typename ValueType , class DeviceType >
class MDArrayDeepCopy< ValueType , DeviceType , DeviceType ,
                       true , true , true >
{
public:

  typedef MDArrayView<ValueType,DeviceType> array_type ;

  static void run( const array_type & dst , const array_type & src )
  {
    typedef Impl::DeepCopyContiguous<ValueType,DeviceType> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

//----------------------------------------------------------------------------
/** \brief  Deep copy from Device to Host,
 *          with different maps,
 *          and no assumption of contiguity.
 *          Force remap to occur on the host.
 */
template< typename ValueType , class DeviceSrc , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceHost , DeviceSrc ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:
  typedef typename DeviceHost::mdarray_map  map_dst ;
  typedef typename DeviceSrc ::mdarray_map  map_src ;
public:

  typedef Serial< HostMemory , map_src >  DeviceTmp ;

  typedef MDArrayView<ValueType,DeviceHost> dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>  tmp_type ;
  typedef MDArrayView<ValueType,DeviceSrc>  src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp ,  DeviceSrc >::run( tmp , src );
    MDArrayDeepCopy< ValueType , DeviceHost , DeviceTmp >::run( dst , tmp );
  }
};

//----------------------------------------------------------------------------
/** \brief  Deep copy from Host to Device,
 *          with different maps,
 *          and no assumption of contiguity.
 *          Force remap to occur on the host.
 */
template< typename ValueType , class DeviceDst , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceDst , DeviceHost ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:
  typedef typename DeviceDst ::mdarray_map  map_dst ;
  typedef typename DeviceHost::mdarray_map  map_src ;
public:

  typedef Serial< HostMemory , map_dst >  DeviceTmp ;

  typedef MDArrayView<ValueType,DeviceDst>  dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>  tmp_type ;
  typedef MDArrayView<ValueType,DeviceHost> src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp , DeviceHost >::run( tmp , src );
    MDArrayDeepCopy< ValueType , DeviceDst , DeviceTmp > ::run( dst , tmp );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MDARRAYVIEW_HPP */


