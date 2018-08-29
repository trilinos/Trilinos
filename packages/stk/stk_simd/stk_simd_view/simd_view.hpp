// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_SIMD_VIEW_H
#define STK_SIMD_VIEW_H

#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_INVALID_INDEX
#define KOKKOS_INVALID_INDEX 0
#endif
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>
#include <stk_simd/Simd.hpp>
#include <stk_simd/Traits.hpp>
#include <stk_simd_view/simd_layout.hpp>
#include <stk_simd_view/simd_index.hpp>
#include <typeinfo>

#define STK_LAMBDA KOKKOS_LAMBDA
#define STK_INLINE KOKKOS_INLINE_FUNCTION
#define STK_FORCE_INLINE KOKKOS_FORCEINLINE_FUNCTION

namespace stk {
namespace simd {

template<typename T>
struct remove_pointer {
  typedef typename std::remove_const<T>::type type;
};

template<typename T>
struct remove_pointer<T*> {
  typedef typename std::remove_const<typename remove_pointer<T>::type>::type type;
};

template<typename T>
struct remove_pointer<T[]> {
  typedef typename std::remove_const<typename remove_pointer<T>::type>::type type;
};

template<typename T, int N>
struct remove_pointer<T[N]> {
  typedef typename std::remove_const<typename remove_pointer<T>::type>::type type;
};

template <class DataType, class ... Properties>
struct ViewTraits;

template <>
struct ViewTraits< void >
{
  typedef void array_layout;
  typedef void execution_space;
  typedef void memory_traits;
};

template< class ... Prop >
struct ViewTraits< void , void , Prop ... >
{
  // Ignore an extraneous 'void'
  typedef typename ViewTraits<void, Prop...>::array_layout    array_layout ;
  typedef typename ViewTraits<void, Prop...>::execution_space execution_space ;
  typedef typename ViewTraits<void, Prop...>::memory_traits   memory_traits ;
};

template< class ArrayLayout , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_array_layout<ArrayLayout>::value >::type , ArrayLayout , Prop ... >
{
  // Specify layout, keep subsequent space and memory traits arguments
  typedef          ArrayLayout                                array_layout ;
  typedef typename ViewTraits<void, Prop...>::execution_space execution_space ;
  typedef typename ViewTraits<void, Prop...>::memory_traits   memory_traits ;
};

template< class Space , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_space<Space>::value >::type , Space , Prop ... >
{
  // Specify Space, memory traits should be the only subsequent argument.
  typedef void                                              array_layout ;
  typedef Space                                             execution_space ;
  typedef typename ViewTraits<void, Prop...>::memory_traits memory_traits ;
};

template< class MemoryTraits , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_memory_traits<MemoryTraits>::value >::type , MemoryTraits , Prop ... >
{
  // Specify memory trait, should not be any subsequent arguments
  typedef void         array_layout ;
  typedef void         execution_space ;
  typedef MemoryTraits memory_traits ;
};

template <class DataType, class ... Prop>
struct ViewTraits {

  typedef ViewTraits< void , Prop ... >  prop ;

  typedef Kokkos::DefaultExecutionSpace default_execution_space;
  typedef Kokkos::MemoryManaged default_memory_traits;
  
  typedef typename remove_pointer<DataType>::type base_type;

  typedef typename
  std::conditional< std::is_same< typename prop::execution_space , void >::value
                    , default_execution_space
                    , typename prop::execution_space
                    >::type
  ExecutionSpace;

    typedef typename
  std::conditional< std::is_same< typename prop::memory_traits , void >::value
                    , default_memory_traits
                    , typename prop::memory_traits
                    >::type
  MemoryTraits;

  // figure out default layout based on execution space
#ifdef KOKKOS_HAVE_CUDA
  static constexpr bool is_device_gpu = std::is_same<ExecutionSpace , Kokkos::Cuda>::value;
#else
  static constexpr bool is_device_gpu = false;
#endif

  typedef typename
  std::conditional< is_device_gpu
                    , Kokkos::LayoutLeft
                    , stk::simd::LayoutRight<base_type> >::type
  default_array_layout;
  
  typedef typename
  std::conditional< std::is_same< typename prop::array_layout , void >::value
                    , default_array_layout
                    , typename prop::array_layout
                    >::type
  ArrayLayout;
};

/*
 * \brief Specialization of Kokkos::View for holding simd types
 *
*/
template< class DataType, class ... Prop >
class View
  : public Kokkos::View< DataType,
                         typename ViewTraits<DataType,Prop...>::ArrayLayout,
                         typename ViewTraits<DataType,Prop...>::ExecutionSpace,
                         typename ViewTraits<DataType,Prop...>::MemoryTraits> {

  typedef stk::simd::View<DataType,Prop...> my_type;

  typedef Kokkos::View< DataType,
                        typename ViewTraits<DataType,Prop...>::ArrayLayout,
                        typename ViewTraits<DataType,Prop...>::ExecutionSpace,
                        typename ViewTraits<DataType,Prop...>::MemoryTraits> view_type;

  typedef typename ViewTraits<DataType,Prop...>::ArrayLayout array_layout;
  typedef typename ViewTraits<DataType,Prop...>::ExecutionSpace execution_space;
  typedef typename ViewTraits<DataType,Prop...>::MemoryTraits memory_traits;

  typedef typename Kokkos::View< DataType, array_layout, execution_space, memory_traits >::traits traits;
  typedef typename Kokkos::View< DataType, array_layout, execution_space, memory_traits >::reference_type reference_type;

  typedef typename stk::SimdT<reference_type>::type simd_reference_type;

  static constexpr bool is_device_gpu = ViewTraits<DataType,Prop...>::is_device_gpu;
  static constexpr int rank = view_type::rank;
  
  typedef typename remove_pointer<DataType>::type base_type;

  static constexpr int simd_width = SimdSizeTraits<base_type>::simd_width;

  static constexpr bool is_layout_right_simd = std::is_same<typename traits::array_layout , stk::simd::LayoutRight<base_type> >::value;
  static constexpr bool is_layout_left_simd = std::is_same<typename traits::array_layout , stk::simd::LayoutLeft<base_type> >::value;
  static constexpr bool is_valid_simd_layout = is_layout_right_simd || is_layout_left_simd || is_device_gpu;

 public:

  using view_type::operator();
  
  STK_INLINE
  my_type & DownCast() const { return static_cast< my_type & > (*this); }

  STK_INLINE
  const my_type & ConstDownCast() const { return static_cast< const my_type & > (*this); }

  //----------------------------------------
  /** \brief  Compatible view of array of scalar types */
  typedef View< typename traits::scalar_array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > 
    array_type ;

  /** \brief  Compatible view of const data type */
  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > 
    const_type ;

  /** \brief  Compatible view of non-const data type */
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > 
    non_const_type ;

  /** \brief  Compatible HostMirror view */
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space >
    HostMirror ;

  STK_INLINE View() : Kokkos::View< DataType, array_layout, execution_space, memory_traits >() {}


  STK_INLINE
  View( const View & rhs ) = default;

  STK_INLINE
  View( View && rhs ) = default;

  STK_INLINE
  View & operator = ( const View & rhs ) = default;

  STK_INLINE
  View & operator = ( View && rhs ) = default;

  //----------------------------------------
  // Compatible view copy constructor and assignment
  // may assign unmanaged from managed.

  template< class RT , class ... RP >
  STK_INLINE
  View( const View<RT,RP...> & rhs )
    : view_type( rhs.ConstDownCast() )
    {}

  template< class RT , class ... RP >
  STK_INLINE
  View & operator = (const View<RT,RP...> & rhs )
  {
    view_type::operator = ( rhs.ConstDownCast() );
    return *this;
  }

  // Allocate with label and layout
  template< typename Label >
  explicit inline
  View( const Label & arg_label,
        typename traits::array_layout const & arg_layout
    )
    : Kokkos::View< DataType, array_layout, execution_space, memory_traits > ( arg_label , arg_layout )
    {}

  template< typename Label >
  explicit STK_INLINE
  View( const Label & arg_label
      , const size_t arg_N0 = KOKKOS_INVALID_INDEX 
      , const size_t arg_N1 = KOKKOS_INVALID_INDEX
      , const size_t arg_N2 = KOKKOS_INVALID_INDEX
      , const size_t arg_N3 = KOKKOS_INVALID_INDEX
      , const size_t arg_N4 = KOKKOS_INVALID_INDEX
      , const size_t arg_N5 = KOKKOS_INVALID_INDEX
      , const size_t arg_N6 = KOKKOS_INVALID_INDEX
      , const size_t arg_N7 = KOKKOS_INVALID_INDEX
      )
    : Kokkos::View< DataType, array_layout, execution_space, memory_traits >( arg_label,
                                                                              arg_N0 , arg_N1 , arg_N2 , arg_N3,
                                                                              arg_N4 , arg_N5 , arg_N6 , arg_N7 )
    {}

  ~View() {}  
  
  #define VALID_SIMD_ERROR_MSG "When calling simd::View with a simd::Index on a platform with a simd length not equal to 1, layouts other that simd::LayoutRight<double/float> and simd::LayoutLeft<double/float> result in undefined behavior."
  #define VALID_INDEX_ERROR_MSG "In operator(Index1, Index2, ...) only Index1 may be of simd::Index type.  All other indices must be an integer type."
  #define VALID_NUM_INDICES_ERROR_MSG "Must call simd::View operator(simd::Index, ...) with the same number of arguments as the rank of the simd::View."

  // one argument

  template <typename I0>
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && is_layout_right_simd, simd_reference_type>::type 
  operator() (const I0 i0) const {
    static_assert(rank==1, VALID_NUM_INDICES_ERROR_MSG);
    return stk::simd::simd_ref_cast(this->data()[simd_width*int_index(i0)]);
  }

  template <typename I0>
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && !is_layout_right_simd, simd_reference_type >::type 
  operator() (const I0 i0) const {
    static_assert(rank==1, VALID_NUM_INDICES_ERROR_MSG);
    static_assert(is_valid_simd_layout, VALID_SIMD_ERROR_MSG);
    return stk::simd::simd_ref_cast((*this)(simd_width*int_index(i0)));
  }

  // two arguments
 
  template< typename I0, typename I1 >
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && is_layout_right_simd, simd_reference_type >::type 
  operator() (const I0 i0, const I1 i1) const {
    static_assert(rank==2, VALID_NUM_INDICES_ERROR_MSG);
    static_assert( !std::is_same<I1, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    return stk::simd::simd_ref_cast(this->data()[simd_width*(i1 + this->extent(1)*int_index(i0))]);
  }

  template< typename I0, typename I1 >
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && !is_layout_right_simd, simd_reference_type >::type 
  operator() (const I0 i0, const I1 i1) const {
    static_assert(rank==2, VALID_NUM_INDICES_ERROR_MSG);
    static_assert( !std::is_same<I1, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    static_assert(is_valid_simd_layout, VALID_SIMD_ERROR_MSG);
    return stk::simd::simd_ref_cast((*this)(simd_width*int_index(i0), i1));
  }

  // three arguments
  
  template< typename I0, typename I1, typename I2 >
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && is_layout_right_simd, simd_reference_type >::type 
  operator() (const I0 i0, const I1 i1, const I2 i2) const {
    static_assert(rank==3, VALID_NUM_INDICES_ERROR_MSG);
    static_assert( !std::is_same<I1, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    static_assert( !std::is_same<I2, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    return stk::simd::simd_ref_cast(this->data()[simd_width*(i2 + this->extent(2)*(
                                                             i1 + this->extent(1)*(
                                                             int_index(i0))))]);
  }

  template< typename I0, typename I1, typename I2 >
  STK_FORCE_INLINE
  typename std::enable_if< std::is_same<I0, stk::simd::Index>::value && !is_device_gpu && !is_layout_right_simd, simd_reference_type >::type 
  operator() (const I0 i0, const I1 i1, const I2 i2) const {
    static_assert(rank==3, VALID_NUM_INDICES_ERROR_MSG);
    static_assert( !std::is_same<I1, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    static_assert( !std::is_same<I2, stk::simd::Index>::value, VALID_INDEX_ERROR_MSG);
    static_assert(is_valid_simd_layout, VALID_SIMD_ERROR_MSG);
    return stk::simd::simd_ref_cast((*this)(simd_width*int_index(i0), i1, i2));
  }

  STK_FORCE_INLINE
  size_t simd_dimension() const {
    return is_device_gpu ? this->extent(0) : simd_pad<base_type>( this->extent(0) ) / simd_width;
  }
};

using Kokkos::create_mirror_view;

template <class DataType, class ... Prop>
typename View<DataType, Prop...>::HostMirror create_mirror_view( View<DataType, Prop...> viewArg) {
  return typename View<DataType, Prop...>::HostMirror(viewArg.label(), viewArg.layout());
}

template <typename A, typename B>
inline void deep_copy(A a, B b) {
  Kokkos::deep_copy(a, b);
}

template <typename A>
inline typename A::HostMirror
copy_from_device(A a) {
  typename A::HostMirror aHost = create_mirror_view(a);
  deep_copy(aHost, a);
  return aHost;
}

}}

#endif // #ifndef SIMD_VIEW__

