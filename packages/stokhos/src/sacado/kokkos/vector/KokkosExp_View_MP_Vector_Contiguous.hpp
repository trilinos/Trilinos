// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_EXPERIMENTAL_VIEW_MP_VECTOR_CONTIGUOUS_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_MP_VECTOR_CONTIGUOUS_HPP

// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_View_MP_Vector_Fwd.hpp"
// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#include "Kokkos_Layout.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

#include "Kokkos_View_Utils.hpp"
#include "Kokkos_View_MP_Vector_Utils.hpp"

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ViewMPVectorContiguous {};

template< class ... Args >
struct is_ViewMPVectorContiguous { enum { value = false }; };

template< class D , class ... P , class ... Args >
struct is_ViewMPVectorContiguous< Kokkos::View<D,P...> , Args... > {
  enum { value =
    std::is_same< typename Kokkos::ViewTraits<D,P...>::specialize
                , ViewMPVectorContiguous >::value
    &&
    ( ( sizeof...(Args) == 0 ) ||
      is_ViewMPVectorContiguous< Args... >::value ) };
};

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

namespace Kokkos {

template <typename T, typename ... P>
struct is_view_mp_vector< View<T,P...> > {
  typedef View<T,P...> view_type;
  static const bool value =
    std::is_same< typename view_type::specialize,
                  Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value;
};

template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr typename
std::enable_if< is_view_mp_vector< View<T,P...> >::value, unsigned >::type
dimension_scalar(const View<T,P...>& view) {
  return view.impl_map().dimension_scalar();
}

} // namespace Kokkos

//----------------------------------------------------------------------------

#include "Sacado_Traits.hpp"
#include "Sacado_MP_Vector.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Kokkos_Core.hpp"

namespace Kokkos {

template <typename D, typename ... P>
struct FlatArrayType< View<D,P...>,
                      typename std::enable_if< is_view_mp_vector< View<D,P...> >::value >::type > {
  typedef View<D,P...> view_type;
  typedef typename view_type::traits::dimension dimension;
  typedef typename view_type::array_type::value_type flat_value_type;
  typedef typename Kokkos::Impl::ViewDataType< flat_value_type , dimension >::type flat_data_type;
  typedef View<flat_data_type,P...> type;
};

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(
  const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  static_assert(std::is_same_v<typename ViewTraits<T, P...>::array_layout, LayoutLeft> ||
    std::is_same_v<typename ViewTraits<T, P...>::array_layout, LayoutRight> ||
    std::is_same_v<typename ViewTraits<T, P...>::array_layout, LayoutStride>);

  using src_type = View<T, P...>;

  auto layout = [&] () {
    if constexpr ( ! std::is_same_v<typename ViewTraits<T, P...>::array_layout, LayoutStride>) {
      return src.layout();
    } else {
      LayoutStride layout;

      for (int idx = 0; idx <= 7; ++idx) {
        layout.dimension[idx] = src.extent(idx);
        layout.stride   [idx] = src.stride(idx);
      }

      return layout;
    }
  }();

  layout.dimension[src_type::rank] = dimension_scalar(src);
  
  const auto prop_copy = Impl::with_properties_if_unset(
    arg_prop, std::string(src.label()).append("_mirror"));

  if constexpr (Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space){
    return typename Impl::MirrorViewType<typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space, T, P ...>::dest_view_type(prop_copy, layout);
  } else {
    return typename View<T, P...>::HostMirror(prop_copy, layout);
  }
}

template <class T, class... P>
inline auto create_mirror(
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  return create_mirror(view_alloc(), src);
}

template <class Space, class T, class... P, typename Enable>
inline auto create_mirror(
  const Space& space,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  return create_mirror(view_alloc(space), src);
}

template <class T, class... P>
inline auto create_mirror(
  Impl::WithoutInitializing_t wi,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  return create_mirror(view_alloc(wi), src);
}

template <class Space, class T, class... P, typename Enable>
inline auto create_mirror(
  Impl::WithoutInitializing_t wi,
  const Space& space,
  const View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  return create_mirror(view_alloc(wi, space), src);
}

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror_view(
  const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
  const Kokkos::View<T, P...>& src,
  typename std::enable_if_t<
    std::is_same_v<typename ViewTraits<T, P...>::specialize,
      Experimental::Impl::ViewMPVectorContiguous>>*
)
{
  return Impl::create_mirror_view(src, arg_prop);
}

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name,
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize,
            Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value &&
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type*)
{
  (void)name;
  fence(
      "Kokkos::create_mirror_view_and_copy: fence before returning src view");  // same behavior as deep_copy(src, src)
  return src;
}

template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name,
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize,
            Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value &&
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type*)
{
  using src_type    = View<T,P...>;
  using Mirror      = typename Impl::MirrorViewType<Space, T, P...>::view_type;
  std::string label = name.empty() ? src.label() : name;
  typename src_type::array_layout layout = src.layout();
  layout.dimension[src_type::rank] = Kokkos::dimension_scalar(src);
  auto mirror       = typename Mirror::non_const_type{
      view_alloc(WithoutInitializing, label), layout};
  deep_copy(mirror, src);
  return mirror;
}

// Overload of deep_copy for MP::Vector views intializing to a constant scalar
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  typedef typename FlatArrayType< View<DT,DP...> >::type flat_array_type;
  Kokkos::Impl::StokhosViewFill< flat_array_type >( view , value );
}

// Overload of deep_copy for MP::Vector views intializing to a constant MP::Vector
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  // static_assert(
  //   Sacado::StaticSize< typename View<DT,DP...>::value_type >::value
  //   ||
  //   std::is_same< Kokkos::Impl::ActiveExecutionMemorySpace
  //               , Kokkos::HostSpace >::value
  //   , "Deep copy from a FAD type must be statically sized or host space" );

  Kokkos::Impl::StokhosViewFill< View<DT,DP...> >( view , value );
}

// Overload of deep_copy for MP::Vector views intializing to a constant scalar
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  typedef typename FlatArrayType< View<DT,DP...> >::type flat_array_type;
  Kokkos::Impl::StokhosViewFill< flat_array_type >( view , value );
}

// Overload of deep_copy for MP::Vector views intializing to a constant MP::Vector
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Can only deep copy into non-const type" );

  // static_assert(
  //   Sacado::StaticSize< typename View<DT,DP...>::value_type >::value
  //   ||
  //   std::is_same< Kokkos::Impl::ActiveExecutionMemorySpace
  //               , Kokkos::HostSpace >::value
  //   , "Deep copy from a FAD type must be statically sized or host space" );

  Kokkos::Impl::StokhosViewFill< View<DT,DP...> >( view , value );
}

/* Specialize for deep copy of MP::Vector */
template< class ExecSpace, class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const ExecSpace & exec,
                const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "Deep copy destination must be non-const" );

  static_assert(
    ( unsigned(ViewTraits<DT,DP...>::rank) ==
      unsigned(ViewTraits<ST,SP...>::rank) )
    , "Deep copy destination and source must have same rank" );

  // Note ETP 09/29/2016:  Use FlatArrayType instead of array_type to work
  // around issue where dst and src are rank-1, but have differing layouts.
  // Kokkos' deep_copy() doesn't work in this case because the array_type
  // will be rank-2.  It should be possible to make deep_copy() work there,
  // but this seems easier.

  // Kokkos::deep_copy(
  //   ExecSpace() ,
  //   typename View<DT,DP...>::array_type( dst ) ,
  //   typename View<ST,SP...>::array_type( src ) );

  Kokkos::deep_copy(
    exec ,
    typename FlatArrayType< View<DT,DP...> >::type( dst ) ,
    typename FlatArrayType< View<ST,SP...> >::type( src ) );
}

/* Specialize for deep copy of MP::Vector */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )>::type * )
{
  using exec_space = typename View<DT,DP...>::execution_space;
  Kokkos::fence();
  Kokkos::deep_copy(exec_space(), dst, src);
  Kokkos::fence();
}

namespace Impl {

template <unsigned N, typename... Args>
KOKKOS_FUNCTION std::enable_if_t<
    N == View<Args...>::rank &&
    std::is_same<typename ViewTraits<Args...>::specialize,
                 Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value,
    View<Args...>>
as_view_of_rank_n(View<Args...> v) {
  return v;
}

// Placeholder implementation to compile generic code for DynRankView; should
// never be called
template <unsigned N, typename T, typename... Args>
std::enable_if_t<
    N != View<T, Args...>::rank &&
        std::is_same<typename ViewTraits<T, Args...>::specialize,
                     Kokkos::Experimental::Impl::ViewMPVectorContiguous>::value,
    View<typename RankDataType<typename View<T, Args...>::value_type, N>::type,
         Args...>>
as_view_of_rank_n(View<T, Args...>) {
  Kokkos::Impl::throw_runtime_exception(
      "Trying to get at a View of the wrong rank");
  return {};
}

}

}

namespace Kokkos {
namespace Impl {

template< class DataType , class ArrayLayout , typename StorageType >
struct ViewDataAnalysis< DataType     /* Original view data type */
                         , ArrayLayout
                         , Sacado::MP::Vector< StorageType > >
{
private:

  typedef typename StorageType::value_type ScalarType;
  typedef ViewArrayAnalysis< DataType > array_analysis ;
  static const int DimVector = StorageType::static_size;

public:

  // Specialized view data mapping:
  typedef Kokkos::Experimental::Impl::ViewMPVectorContiguous specialize ;

  typedef typename array_analysis::dimension             dimension ;
  typedef typename array_analysis::value_type            value_type ;
  typedef typename array_analysis::const_value_type      const_value_type ;
  typedef typename array_analysis::non_const_value_type  non_const_value_type ;

  // Generate analogous multidimensional array specification type.
  typedef typename
    ViewDataType< value_type , dimension >::type  type ;
  typedef typename
    ViewDataType< const_value_type , dimension >::type  const_type ;
  typedef typename
    ViewDataType< non_const_value_type , dimension >::type  non_const_type ;

private:

  // A const ?
  enum { is_const = std::is_same< value_type , const_value_type >::value };

  // The unwrapped scalar types:
  typedef typename
    std::conditional< is_const , const ScalarType , ScalarType >::type
      scalar_type ;

  typedef ScalarType        non_const_scalar_type ;
  typedef const ScalarType  const_scalar_type ;

  // Prepend or append the vector dimension based on ArrayLayout
  // Note:  you can't prepend a static dimension, so use 0 for LayoutLeft
  typedef typename array_analysis::dimension::
    template prepend<0>::type
      prepend_scalar_dimension ;
  typedef typename array_analysis::dimension::
    template append<DimVector>::type
      append_scalar_dimension ;
  typedef typename std::conditional<
    std::is_same< ArrayLayout, Kokkos::LayoutLeft>::value,
    prepend_scalar_dimension,
    append_scalar_dimension >::type scalar_dimension;

public:

  // Generate "flattened" multidimensional array specification type.
  typedef typename
    ViewDataType< scalar_type , scalar_dimension >::type scalar_array_type ;

  typedef typename
    ViewDataType< const_scalar_type , scalar_dimension >::type
      const_scalar_array_type ;

  typedef typename
    ViewDataType< non_const_scalar_type , scalar_dimension >::type
      non_const_scalar_array_type ;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

  template < class ValueType,
             bool is_static = Sacado::IsStaticallySized<ValueType>::value >
struct MPVectorAllocation;

// MP::Vector allocation for statically-sized MP::Vector types.
// In this case we can reinterpret cast directly between pointers of types
// MP::Vector<Storage> and MP::Vector<Storage>::value_type.
template <class ValueType>
struct MPVectorAllocation<ValueType, true> {
  typedef ValueType value_type;
  typedef typename Sacado::ValueType<value_type>::type scalar_type;

  value_type  * value_ptr;
  scalar_type * scalar_ptr;

  KOKKOS_INLINE_FUNCTION
  static constexpr size_t
  memory_span(const size_t span, const unsigned vector_size) {
    return span * vector_size * sizeof(scalar_type);
  }

  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation() : value_ptr(0), scalar_ptr(0) {}

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation& operator=(const MPVectorAllocation<T,true>& a) {
    value_ptr = a.value_ptr;
    scalar_ptr = a.scalar_ptr;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void set(value_type* ptr, const size_t span, const unsigned vector_size) {
    value_ptr  = ptr;
    scalar_ptr = reinterpret_cast<scalar_type*>(ptr);
  }

  template <class ExecSpace>
  struct ConstructDestructFunctor {
    typedef Kokkos::Impl::ViewValueFunctor< ExecSpace, scalar_type > FunctorType ;
    FunctorType m_functor;
    bool m_initialize;

    ConstructDestructFunctor() = default;
    ConstructDestructFunctor(const ConstructDestructFunctor&) = default;
    ConstructDestructFunctor& operator=(const ConstructDestructFunctor&) = default;

    ConstructDestructFunctor(const ExecSpace & space,
                             const bool initialize,
                             const size_t span,
                             const unsigned vector_size,
                             scalar_type* scalar_ptr) :
      m_functor( space , scalar_ptr , span*vector_size , "Stokhos_MP_VectorContig_ConstructDestructFunctor1" ),
      m_initialize(initialize) {}

    inline void construct_shared_allocation() {
      if (m_initialize)
        m_functor.construct_shared_allocation();
    }

    inline void destroy_shared_allocation() {
      if (m_initialize)
        m_functor.destroy_shared_allocation();
    }

  };

  template <class ExecSpace>
  inline ConstructDestructFunctor<ExecSpace>
  create_functor(const ExecSpace & space,
                 const bool initialize,
                 const size_t span,
                 const unsigned vector_size) const {
    return ConstructDestructFunctor<ExecSpace>(space, initialize, span, vector_size, scalar_ptr);
  }

  // Assign scalar_type pointer to give ptr
  template <typename T>
  KOKKOS_INLINE_FUNCTION
  void assign(T * ptr) {
    value_ptr  = reinterpret_cast<value_type*>(ptr);
    scalar_ptr = reinterpret_cast<scalar_type*>(ptr);
  }

};

// MP::Vector allocation for dynamically-sized MP::Vector types.
// In this case we allocate two chunks of data, the first for the the
// MP::Vector<Storage> itself and then for the underlying scalar type
// (MP::Vector<Storage>::value_type).  The memory is laid out with the
// former followed by the latter.
template <class ValueType>
struct MPVectorAllocation<ValueType, false> {
  typedef ValueType value_type;
  typedef typename Sacado::ValueType<value_type>::type scalar_type;

  value_type  * value_ptr;
  scalar_type * scalar_ptr;

  KOKKOS_INLINE_FUNCTION
  static constexpr size_t
  memory_span(const size_t span, const unsigned vector_size) {
    return span * ( vector_size * sizeof(scalar_type) + sizeof(value_type) );
  }

  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation() : value_ptr(0), scalar_ptr(0) {}

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  MPVectorAllocation& operator=(const MPVectorAllocation<T,false>& a) {
    value_ptr = a.value_ptr;
    scalar_ptr = a.scalar_ptr;
    return *this;
  }

  // We are making an assumption the data is laid out as described above,
  // which in general may not be true if the view is created from memory
  // allocated elsewhere.  We should check for that.
  KOKKOS_INLINE_FUNCTION
  void set(value_type* ptr, const size_t span, const unsigned vector_size) {
    value_ptr = ptr;
    scalar_ptr = reinterpret_cast<scalar_type*>(ptr+span);
  }

  template <class ExecSpace>
  struct VectorConstruct {
    ExecSpace m_space;
    value_type* m_p;
    scalar_type* m_sp;
    size_t m_span;
    unsigned m_vector_size;

    VectorConstruct() = default;
    VectorConstruct(const VectorConstruct&) = default;
    VectorConstruct& operator=(const VectorConstruct&) = default;

    inline
    VectorConstruct(const ExecSpace& space,
                    value_type* p,
                    scalar_type* sp,
                    const size_t span,
                    const unsigned vector_size) :
      m_space(space), m_p(p), m_sp(sp), m_span(span), m_vector_size(vector_size) {}

    inline void execute() {
      typedef Kokkos::RangePolicy< ExecSpace > PolicyType ;
      const Kokkos::Impl::ParallelFor< VectorConstruct , PolicyType >
        closure( *this , PolicyType( 0 , m_span ) );
      closure.execute();
      m_space.fence();
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_t i) const {
      new (m_p+i) value_type(m_vector_size, m_sp+i*m_vector_size, false);
    }
  };

  template <class ExecSpace>
  struct ConstructDestructFunctor {
    typedef Kokkos::Impl::ViewValueFunctor< ExecSpace, scalar_type > ScalarFunctorType ;
    typedef VectorConstruct< ExecSpace > VectorFunctorType ;
    ScalarFunctorType m_scalar_functor;
    VectorFunctorType m_vector_functor;
    bool m_initialize;

    ConstructDestructFunctor() = default;
    ConstructDestructFunctor(const ConstructDestructFunctor&) = default;
    ConstructDestructFunctor& operator=(const ConstructDestructFunctor&) = default;

    ConstructDestructFunctor(const ExecSpace & space,
                             const bool initialize,
                             const size_t span,
                             const unsigned vector_size,
                             scalar_type* scalar_ptr,
                             value_type* value_ptr) :
      m_scalar_functor( space , scalar_ptr , span*vector_size , "Stokhos_MP_VectorContig_ConstructDestructFunctor2" ),
      m_vector_functor( space , value_ptr , scalar_ptr , span , vector_size ),
      m_initialize(initialize) {}

    inline void construct_shared_allocation() {
      // First initialize the scalar_type array
      if (m_initialize)
        m_scalar_functor.construct_shared_allocation();

      // Construct each MP::Vector using memory in scalar_ptr array,
      // setting pointer to MP::Vector values from values array
      // Equivalent to:
      // value_type* p = value_ptr;
      // scalar_type* sp = scalar_ptr;
      // for (size_t i=0; i<span; ++i) {
      //   new (p++) value_type(vector_size, sp, false);
      //   sp += vector_size;
      // }
      // (we always need to do this, regardless of initialization)
      m_vector_functor.execute();
    }

    inline void destroy_shared_allocation() {
      // We only need to (possibly) call the destructor on values in the
      // scalar_type array, since the value_type array is a view into it
      if (m_initialize)
        m_scalar_functor.destroy_shared_allocation();
    }

  };

  template <class ExecSpace>
  inline ConstructDestructFunctor<ExecSpace>
  create_functor(const ExecSpace & space,
                 const bool initialize,
                 const size_t span,
                 const unsigned vector_size) const {
    return ConstructDestructFunctor<ExecSpace>(space, initialize, span, vector_size, scalar_ptr, value_ptr);
  }

  // Assign scalar_type pointer to give ptr
  // This makes BIG assumption on how the data was allocated
  template <typename T>
  KOKKOS_INLINE_FUNCTION
  void assign(T * ptr) {
    value_ptr  = reinterpret_cast<value_type*>(ptr);
    if (ptr != 0)
      scalar_ptr = value_ptr->coeff();
    else
      scalar_ptr = 0;
  }
};

}}} // namespace Kokkos::Experimental::Impl

namespace Kokkos {
namespace Impl {

template< class Traits >
class ViewMapping< Traits , /* View internal mapping */
  typename std::enable_if<
    ( std::is_same< typename Traits::specialize
                  , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
      &&
      ( std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutLeft >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutRight >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    )
    , typename Traits::specialize
    >::type >
{
private:

  template< class , class ... > friend class ViewMapping ;
  template< class , class ... > friend class Kokkos::View ;

  typedef typename Traits::value_type  sacado_mp_vector_type ;
  typedef typename sacado_mp_vector_type::storage_type stokhos_storage_type ;
  typedef typename stokhos_storage_type::value_type intrinsic_scalar_type ;
  typedef typename
    std::add_const< intrinsic_scalar_type >::type  const_intrinsic_scalar_type ;

  enum { StokhosStorageStaticDimension = stokhos_storage_type::static_size };
  typedef Sacado::integral_nonzero< unsigned , StokhosStorageStaticDimension > sacado_size_type;

  typedef Kokkos::Experimental::Impl::MPVectorAllocation<sacado_mp_vector_type> handle_type;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  // Prepend or append the vector dimension based on array_layout
  // Note:  you can't prepend a static dimension, so use 0 for LayoutLeft
  typedef ViewArrayAnalysis< typename Traits::data_type > array_analysis ;
  typedef typename array_analysis::dimension array_dimension;
  typedef ViewOffset< typename array_dimension::
                        template append<StokhosStorageStaticDimension>::type,
                      typename Traits::array_layout,
                      void
                      >  append_offset_type ;
  typedef ViewOffset< typename array_dimension::
                        template prepend<0>::type,
                      typename Traits::array_layout,
                      void
                      >  prepend_offset_type ;
  typedef typename std::conditional<
    std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft>::value,
    prepend_offset_type,
    append_offset_type >::type array_offset_type;

  handle_type      m_impl_handle ;
  offset_type      m_impl_offset ;
  unsigned         m_stride ;
  sacado_size_type m_sacado_size ; // Size of sacado dimension

  // Note:  if the view is partitioned, m_sacado_size is not the stride in
  // memory between consecutive MP::Vector entries for given vector index:
  //
  // original_sacado_size = m_stride * m_sacado_size
  // m_stride = 1 for original allocation.
  //
  // Stride here has a slightly different meaning than in the standard
  // View implementation.  For the moment we are assuming no padding within
  // the view array itself and stride is to allow for partitioning the view
  // by dividing up the scalar type.
  //
  // I suspect we could combine this with the way the stride is managed in
  // the default view, in which case, I don't think we even need a
  // specialization
  //
  // For reshaping by folding the sacado dimension into its next adjacent
  // dimension, padding wouldn't generally work.  So unless there becomes
  // a way to turn padding off in the default view, a specialization
  // will be necessary.

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  // Rank corresponding to the sacado dimension
  enum { Sacado_Rank = std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft >::value ? 0 : Rank+1 };

  // Using the internal offset mapping so limit to public rank:
  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr size_t extent( const iType & r ) const
    { return m_impl_offset.m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr
  typename Traits::array_layout layout() const
    { return m_impl_offset.layout(); }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const
    { return m_impl_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const
    { return m_impl_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const
    { return m_impl_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const
    { return m_impl_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const
    { return m_impl_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const
    { return m_impl_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const
    { return m_impl_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const
    { return m_impl_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  // Since we all for striding within the data type, we can't guarantee
  // regular striding
  using is_regular = std::false_type ;

  // FIXME:  Adjust these for m_stride
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const
    { return m_impl_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const
    { return m_impl_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const
    { return m_impl_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const
    { return m_impl_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const
    { return m_impl_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const
    { return m_impl_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const
    { return m_impl_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const
    { return m_impl_offset.stride_7(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const
    { m_impl_offset.stride(s); }

  // Size of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned dimension_scalar() const
    { return m_sacado_size.value; }

  // Whether the storage type is statically sized
  static const bool is_static = stokhos_storage_type::is_static ;

  // Whether sacado dimension is contiguous
  static const bool is_contiguous = true;

  //----------------------------------------
  // Range of mapping

  // Return type of reference operators
  typedef sacado_mp_vector_type & reference_type ;

  /** \brief Pointer to underlying memory type */
  typedef sacado_mp_vector_type * pointer_type ;

  /** \brief  Span of the mapped range : [ data() .. data() + span() ) */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const
    { return m_impl_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_impl_offset.span_is_contiguous() && (m_stride == 1); }

  /** \brief Raw data access */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const
    { return m_impl_handle.value_ptr ; }

  //----------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const
    { return *m_impl_handle.value_ptr; }

  // FIXME:  Check this
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    ! std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const
    { return m_impl_handle.value_ptr[m_stride * i0]; }

  // FIXME:  Check this
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0) ]; }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1) ]; }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2,i3) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2,i3,i4) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2,i3,i4,i5) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2,i3,i4,i5,i6) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 ) const
    { return m_impl_handle.value_ptr[ m_stride * m_impl_offset(i0,i1,i2,i3,i4,i5,i6,i7) ]; }

  //----------------------------------------

  /** \brief  Span, in bytes, of the required memory */
  KOKKOS_INLINE_FUNCTION
  static size_t memory_span( typename Traits::array_layout const & layout )
    {
      // Do not introduce padding...
      typedef std::integral_constant< unsigned , 0 >  padding ;
      offset_type offset( padding(), layout );

      // Always use static dimension if we are static
      const unsigned static_dim = StokhosStorageStaticDimension;
      if (static_dim > 0)
        return handle_type::memory_span( offset.span(), static_dim );

      // Else get size from prescribed layout
      const size_t sacado_size =
        Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout);
      return handle_type::memory_span( offset.span(), sacado_size );
    }

  //----------------------------------------

  KOKKOS_DEFAULTED_FUNCTION ~ViewMapping() = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping() :
    m_impl_handle(),
    m_impl_offset(),
    m_stride(1),
    m_sacado_size(0)
    {}

  KOKKOS_DEFAULTED_FUNCTION ViewMapping( const ViewMapping & ) = default ;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping & operator = ( const ViewMapping & ) = default ;

  KOKKOS_DEFAULTED_FUNCTION ViewMapping( ViewMapping && ) = default ;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping & operator = ( ViewMapping && ) = default ;

  template< class ... P >
  KOKKOS_INLINE_FUNCTION
  ViewMapping
    ( ViewCtorProp< P ... > const & prop
    , typename Traits::array_layout const & layout
    )
    : m_impl_handle()
    , m_impl_offset( std::integral_constant< unsigned , 0 >()
              , layout )
    , m_stride( 1 )
    , m_sacado_size( Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout) )
    {
      m_impl_handle.set( ( (ViewCtorProp<void,pointer_type> const &) prop ).value,
                    m_impl_offset.span(), m_sacado_size.value );
    }

  /**\brief  Assign data */
  KOKKOS_INLINE_FUNCTION
  void assign_data( pointer_type arg_ptr )
  { m_impl_handle.set( arg_ptr, m_impl_offset.span(), m_sacado_size.value ); }

  //----------------------------------------
  /*  Allocate and construct mapped array.
   *  Allocate via shared allocation record and
   *  return that record for allocation tracking.
   */
  template< class ... P >
  SharedAllocationRecord<> *
  allocate_shared( ViewCtorProp< P... > const & prop
                 , typename Traits::array_layout const & layout )
  {
    typedef ViewCtorProp< P... > ctor_prop ;

    typedef typename ctor_prop::execution_space  execution_space ;
    typedef typename Traits::memory_space         memory_space ;
    typedef typename handle_type::template ConstructDestructFunctor<execution_space> functor_type ;
    typedef SharedAllocationRecord< memory_space , functor_type > record_type ;

    // Disallow padding
    typedef std::integral_constant< unsigned , 0 > padding ;

    m_impl_offset = offset_type( padding(), layout );
    m_stride = 1;
    m_sacado_size = Kokkos::Impl::GetSacadoSize<unsigned(Rank)>::eval(layout);

    const size_t alloc_size =
      handle_type::memory_span( m_impl_offset.span(), m_sacado_size.value );

    // Create shared memory tracking record with allocate memory from the memory space
    record_type * const record =
      record_type::allocate( ( (ViewCtorProp<void,memory_space> const &) prop ).value
                           , ( (ViewCtorProp<void,std::string>  const &) prop ).value
                           , alloc_size );

    //  Only set the the pointer and initialize if the allocation is non-zero.
    //  May be zero if one of the dimensions is zero.
    if ( alloc_size ) {

      auto space = ((ViewCtorProp<void,execution_space> const &) prop).value;
      m_impl_handle.set( reinterpret_cast< pointer_type >( record->data() ),
                    m_impl_offset.span(), m_sacado_size.value );

      // Assume destruction is only required when construction is requested.
      // The ViewValueFunctor has both value construction and destruction operators.
      record->m_destroy = m_impl_handle.create_functor(
        space
        , ctor_prop::initialize
        , m_impl_offset.span()
        , m_sacado_size.value );

      // Construct values
      record->m_destroy.construct_shared_allocation();
      space.fence();
    }

    return record ;
  }

  template< class ... P >
  SharedAllocationRecord<> *
  allocate_shared( ViewCtorProp< P... > const & prop
                 , typename Traits::array_layout const & layout
                 , bool /*execution_space_specified*/)
  {
    return allocate_shared(prop, layout);
  }

  //----------------------------------------
  // If the View is to construct or destroy the elements.

  /*
  template< class ExecSpace >
  void construct( const ExecSpace & space ) const
    {
      m_impl_handle.construct( space, m_impl_offset.span(), m_sacado_size.value );
    }

  template< class ExecSpace >
  void destroy( const ExecSpace & space ) const
    {
      m_impl_handle.destruct( space, m_impl_offset.span(), m_sacado_size.value );
    }
  */
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Assign compatible Sacado::MP::Vector view mappings.
 *
 *  View<MP::Vector> = View<MP::Vector>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has MP::Vector
    std::is_same< typename DstTraits::specialize
                , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    &&
    // Source view has MP::Vector only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of MP::Vector requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ||
        ( unsigned(DstTraits::rank) == 0 && unsigned(SrcTraits::rank) == 0 ) ||
        ( unsigned(DstTraits::rank) == 1 && unsigned(SrcTraits::rank) == 1 ) ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable
          < typename DstType::offset_type::dimension_type
          , typename SrcType::offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

      dst.m_impl_handle  = src.m_impl_handle ;
      dst.m_impl_offset  = src.m_impl_offset ;
      dst.m_stride  = src.m_stride ;
      dst.m_sacado_size = src.m_sacado_size ;
    }
};

/**\brief  Assign compatible Sacado::MP::Vector view mappings.
 *
 *  View<ordinary> = View<MP::Vector>
 *  where View<ordinay>::Rank = View<MP::Vector>::Rank+1
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has MP::Vector only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    &&
    // Ranks match
    unsigned(DstTraits::dimension::rank) == unsigned(SrcTraits::dimension::rank)+1
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of MP::Vector requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::scalar_array_type
                    , typename SrcTraits::scalar_array_type >::value ||
        std::is_same< typename DstTraits::scalar_array_type
                    , typename SrcTraits::const_scalar_array_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable<
          typename DstType::offset_type::dimension_type,
          typename SrcType::array_offset_type::dimension_type >::value,
        "View assignment must have compatible dimensions" );

      if ( src.m_stride != 1 ) {
        Kokkos::abort("\n\n ****** Kokkos::View< Sacado::MP::Vector ... > cannot assign with non-unit stride ******\n\n");
      }

      unsigned dims[8];
      dims[0] = src.m_impl_offset.dimension_0();
      dims[1] = src.m_impl_offset.dimension_1();
      dims[2] = src.m_impl_offset.dimension_2();
      dims[3] = src.m_impl_offset.dimension_3();
      dims[4] = src.m_impl_offset.dimension_4();
      dims[5] = src.m_impl_offset.dimension_5();
      dims[6] = src.m_impl_offset.dimension_6();
      dims[7] = src.m_impl_offset.dimension_7();
      unsigned rank = SrcTraits::dimension::rank;
      unsigned sacado_size = src.m_sacado_size.value;
      if (std::is_same<typename SrcTraits::array_layout, LayoutLeft>::value) {
        // Move sacado_size to the first dimension, shift all others up one
        for (unsigned i=rank; i>0; --i)
          dims[i] = dims[i-1];
        dims[0] = sacado_size;
      }
      else {
        dims[rank] = sacado_size;
      }
      typedef typename DstType::offset_type dst_offset_type;
      dst.m_impl_offset = dst_offset_type( std::integral_constant< unsigned , 0 >(),
                                      typename DstTraits::array_layout(
                                        dims[0] , dims[1] , dims[2] , dims[3] ,
                                        dims[4] , dims[5] , dims[6] , dims[7] ) );
      dst.m_impl_handle  = src.m_impl_handle.scalar_ptr ;
    }
};

/**\brief  Assign compatible Sacado::MP::Vector view mappings.
 *
 *  View<ordinary> = View<MP::Vector>
 *  where View<ordinay>::Rank = View<MP::Vector>::Rank, i.e., assigning
 *  to the "flattened" view type
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has MP::Vector only
    std::is_same< typename SrcTraits::specialize
                , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    &&
    // Ranks match
    unsigned(DstTraits::dimension::rank) == unsigned(SrcTraits::dimension::rank)
    )
    , typename DstTraits::specialize
    >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcType ;

  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of MP::Vector requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::non_const_value_type::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , const typename SrcTraits::non_const_value_type::value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable<
          typename DstType::offset_type::dimension_type,
          typename SrcType::offset_type::dimension_type >::value,
        "View assignment must have compatible dimensions" );

      if ( src.m_stride != 1 ) {
       Kokkos::abort("\n\n ****** Kokkos::View< Sacado::MP::Vector ... > cannot assign with non-unit stride ******\n\n");
      }

      unsigned dims[8];
      dims[0] = src.m_impl_offset.dimension_0();
      dims[1] = src.m_impl_offset.dimension_1();
      dims[2] = src.m_impl_offset.dimension_2();
      dims[3] = src.m_impl_offset.dimension_3();
      dims[4] = src.m_impl_offset.dimension_4();
      dims[5] = src.m_impl_offset.dimension_5();
      dims[6] = src.m_impl_offset.dimension_6();
      dims[7] = src.m_impl_offset.dimension_7();
      unsigned rank = SrcTraits::dimension::rank;
      unsigned sacado_size = src.m_sacado_size.value;
      if (std::is_same<typename DstTraits::array_layout, LayoutLeft>::value) {
        dims[0] = dims[0]*sacado_size;
        dims[rank] = 0;
      }
      else {
        dims[rank-1] = dims[rank-1]*sacado_size;
        dims[rank] = 0;
      }
      typedef typename DstType::offset_type dst_offset_type;
      dst.m_impl_offset = dst_offset_type( std::integral_constant< unsigned , 0 >(),
                                      typename DstTraits::array_layout(
                                        dims[0] , dims[1] , dims[2] , dims[3] ,
                                        dims[4] , dims[5] , dims[6] , dims[7] ) );
      dst.m_impl_handle  = src.m_impl_handle.scalar_ptr ;
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Subview mapping

template< class DataType, class ... P , class Arg0, class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      // Source view has MP::Vector only
      std::is_same< typename Kokkos::ViewTraits<DataType,P...>::specialize
                  , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
      &&
      (
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename Kokkos::ViewTraits<DataType,P...>::array_layout
                    , Kokkos::LayoutStride >::value
      )
      && !Sacado::MP::is_vector_partition<Arg0>::value
    )>::type
  , Kokkos::ViewTraits<DataType,P...>
  , Arg0, Args ... >
{
private:

  typedef Kokkos::ViewTraits<DataType,P...> SrcTraits;

  //static_assert( SrcTraits::rank == sizeof...(Args) , "" );

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Arg0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Arg0,Args...>::value)
    , R2 = bool(is_integral_extent<2,Arg0,Args...>::value)
    , R3 = bool(is_integral_extent<3,Arg0,Args...>::value)
    , R4 = bool(is_integral_extent<4,Arg0,Args...>::value)
    , R5 = bool(is_integral_extent<5,Arg0,Args...>::value)
    , R6 = bool(is_integral_extent<6,Arg0,Args...>::value)
    };

  // Public rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  // Whether right-most non-MP::Vector rank is a range.
  enum { R0_rev = ( 0 == SrcTraits::rank ? RZ : (
                    1 == SrcTraits::rank ? R0 : (
                    2 == SrcTraits::rank ? R1 : (
                    3 == SrcTraits::rank ? R2 : (
                    4 == SrcTraits::rank ? R3 : (
                    5 == SrcTraits::rank ? R4 : (
                    6 == SrcTraits::rank ? R5 : R6 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1 or 2, InputLayout Left, Interval 0
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0 && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0_rev && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value )
      ), typename SrcTraits::array_layout , Kokkos::LayoutStride
      >::type array_layout ;

  typedef typename SrcTraits::value_type  sacado_mp_vector_type ;

  typedef typename std::conditional< rank == 0 , sacado_mp_vector_type ,
          typename std::conditional< rank == 1 , sacado_mp_vector_type * ,
          typename std::conditional< rank == 2 , sacado_mp_vector_type ** ,
          typename std::conditional< rank == 3 , sacado_mp_vector_type *** ,
          typename std::conditional< rank == 4 , sacado_mp_vector_type **** ,
          typename std::conditional< rank == 5 , sacado_mp_vector_type ***** ,
          typename std::conditional< rank == 6 , sacado_mp_vector_type ****** ,
                                                 sacado_mp_vector_type *******
          >::type >::type >::type >::type >::type >::type >::type
    data_type ;

public:

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  // The presumed type is 'ViewMapping< traits_type , void >'
  // However, a compatible ViewMapping is acceptable.
  template< class DstTraits >
  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< DstTraits , typename DstTraits::specialize > & dst
                    , ViewMapping< SrcTraits , typename SrcTraits::specialize > const & src
                    , Arg0 arg0, Args ... args )
    {
      static_assert(
        ViewMapping< DstTraits , traits_type , typename DstTraits::specialize >::is_assignable ,
        "Subview destination type must be compatible with subview derived type" );

      typedef ViewMapping< DstTraits , typename DstTraits::specialize > DstType ;
      typedef typename DstType::offset_type  dst_offset_type ;

      const SubviewExtents< SrcTraits::rank , rank >
        extents( src.m_impl_offset.m_dim , arg0 , args... );

      const size_t offset = src.m_impl_offset( extents.domain_offset(0)
                                          , extents.domain_offset(1)
                                          , extents.domain_offset(2)
                                          , extents.domain_offset(3)
                                          , extents.domain_offset(4)
                                          , extents.domain_offset(5)
                                          , extents.domain_offset(6)
                                          , extents.domain_offset(7) );

      dst.m_impl_offset = dst_offset_type( src.m_impl_offset , extents );
      dst.m_impl_handle.value_ptr = src.m_impl_handle.value_ptr + offset;
      dst.m_impl_handle.scalar_ptr =
        src.m_impl_handle.scalar_ptr + offset * src.m_stride * src.m_sacado_size.value;
      dst.m_stride = src.m_stride;
      dst.m_sacado_size = src.m_sacado_size;
    }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Partition mapping

template< class DataType, class ...P, unsigned Size >
class ViewMapping<
  void,
  ViewTraits<DataType,P...> ,
  Sacado::MP::VectorPartition<Size> >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef ViewTraits<DataType,P...> src_traits;
  typedef ViewMapping< src_traits , typename src_traits::specialize >  src_type ;

  typedef typename src_type::offset_type::dimension_type src_dimension;
  typedef typename src_traits::value_type mp_vector_type;
  typedef typename mp_vector_type::storage_type storage_type;
  typedef typename storage_type::template apply_N<Size> storage_apply;
  typedef typename storage_apply::type strided_storage_type;
  typedef Sacado::MP::Vector< strided_storage_type > strided_value_type;
  typedef typename
    ViewDataType< strided_value_type , src_dimension >::type strided_data_type;
  typedef ViewTraits<strided_data_type,P...> dst_traits;
  typedef View<strided_data_type,P...> type;
  typedef ViewMapping< dst_traits , typename dst_traits::specialize >  dst_type ;

  KOKKOS_INLINE_FUNCTION static
  void assign( dst_type & dst
             , const src_type & src
             , const Sacado::MP::VectorPartition<Size> & part )
    {
      // The pointer assignments below are not sufficient for dynamically sized
      // scalar types, so disallow this case for now
      static_assert( storage_type::is_static,
                     "For performance reasons, partitioned assignment is only implemented for statically-sized MP::Vector types" );

      unsigned len = part.end - part.begin;
      if ( Size != len || Size == 0 ) {
        Kokkos::abort("\n\n ******  Kokkos::View< Sacado::MP::Vector ... > Invalid size in partitioned view assignment ******\n\n");
      }

      dst.m_impl_handle.value_ptr =
        reinterpret_cast<strided_value_type*>( src.m_impl_handle.value_ptr ) +
        part.begin / len ;
      dst.m_impl_handle.scalar_ptr = src.m_impl_handle.scalar_ptr +
        (part.begin / len) * src.m_stride * src.m_sacado_size.value ;
      dst.m_impl_offset  = src.m_impl_offset ;
      dst.m_stride  = src.m_stride * src.m_sacado_size.value / Size ;
      dst.m_sacado_size = len ;
    }
};

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template< unsigned Size, typename D, typename ... P  >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::ViewMapping< void, typename Kokkos::ViewTraits<D,P...>, Sacado::MP::VectorPartition<Size> >::type
partition( const Kokkos::View<D,P...> & src ,
           const unsigned beg )
{
  typedef Kokkos::ViewTraits<D,P...> traits;
  typedef typename Kokkos::Impl::ViewMapping< void, traits, Sacado::MP::VectorPartition<Size> >::type DstViewType;
  const Sacado::MP::VectorPartition<Size> part( beg , beg+Size );
  return DstViewType(src, part);
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Specialization for deep_copy( view, view::value_type ) for Cuda
#if defined( KOKKOS_ENABLE_CUDA )
template< class OutputView >
struct StokhosViewFill< OutputView ,
                 typename std::enable_if< std::is_same< typename OutputView::specialize,
                                                        Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value &&
                                     std::is_same< typename OutputView::execution_space,
                                                   Cuda >::value >::type >
{
  typedef typename OutputView::const_value_type   const_value_type ;
  typedef typename OutputView::execution_space    execution_space ;
  typedef typename OutputView::size_type          size_type ;

  template <unsigned VectorLength>
  struct Kernel {
    typedef typename OutputView::execution_space execution_space ;
    const OutputView output;
    const_value_type input;

    Kernel( const OutputView & arg_out , const_value_type & arg_in ) :
      output(arg_out), input(arg_in) {}

    typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member ;

    KOKKOS_INLINE_FUNCTION
    void operator()( const team_member & dev ) const
    {
      const size_type tidx = dev.team_rank() % VectorLength;
      const size_type tidy = dev.team_rank() / VectorLength;
      const size_type nrow = dev.team_size() / VectorLength;
      const size_type nvec = dimension_scalar(output);

      const size_type i0 = dev.league_rank() * nrow + tidy;
      if ( i0 >= output.extent(0) ) return;

      for ( size_type i1 = 0 ; i1 < output.extent(1) ; ++i1 ) {
      for ( size_type i2 = 0 ; i2 < output.extent(2) ; ++i2 ) {
      for ( size_type i3 = 0 ; i3 < output.extent(3) ; ++i3 ) {
      for ( size_type i4 = 0 ; i4 < output.extent(4) ; ++i4 ) {
      for ( size_type i5 = 0 ; i5 < output.extent(5) ; ++i5 ) {
      for ( size_type i6 = 0 ; i6 < output.extent(6) ; ++i6 ) {
      for ( size_type i7 = 0 ; i7 < output.extent(7) ; ++i7 ) {
      for ( size_type is = tidx ; is < nvec ; is+=VectorLength ) {
        output.access(i0,i1,i2,i3,i4,i5,i6,i7).fastAccessCoeff(is) =
          input.fastAccessCoeff(is) ;
      }}}}}}}}
    }
  };

  StokhosViewFill( const OutputView & output , const_value_type & input )
  {
    if ( Sacado::is_constant(input) ) {
      deep_copy( output , input.fastAccessCoeff(0) );
    }
    else {

      // Coalesced accesses are 128 bytes in size
      typedef typename OutputView::array_type::value_type scalar_type;
      const unsigned vector_length =
        ( 128 + sizeof(scalar_type)-1 ) / sizeof(scalar_type);

      // 8 warps per block should give good occupancy
      const size_type block_size = 256;

      const size_type rows_per_block = block_size / vector_length;
      const size_type n = output.extent(0);
      const size_type league_size = ( n + rows_per_block-1 ) / rows_per_block;
      const size_type team_size = rows_per_block * vector_length;
      Kokkos::TeamPolicy< execution_space > config( league_size, team_size );

      parallel_for( config, Kernel<vector_length>(output, input) );
      execution_space().fence();
    }
  }

};
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct ViewSpecializeSacadoFad;

/**\brief  Assign compatible Sacado FAD<MP::Vector> view mappings.
 *
 *  View<MP::Vector> = View< FAD<MP::Vector> >
 *
 * This only works for statically allocated MP::Vector currently
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has MP::Vector only
    std::is_same< typename DstTraits::specialize
                , Kokkos::Experimental::Impl::ViewMPVectorContiguous >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFad >::value
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template< class DstType >
  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcFadType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::data_type
                    , typename SrcTraits::scalar_array_type >::value ||
        std::is_same< typename DstTraits::data_type
                    , typename SrcTraits::const_scalar_array_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable
          < typename DstType::offset_type::dimension_type
          , typename SrcFadType::array_offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

      typedef typename DstType::offset_type  dst_offset_type ;

      dst.m_impl_offset  = dst_offset_type( src.m_array_offset );
      dst.m_impl_handle.assign(src.m_impl_handle) ;
      dst.m_stride  = 1;

      // Don't need to set dst.m_sacado_size since it is determined statically
      static_assert( DstType::is_static,
                     "Destination view must be statically allocated" );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "Kokkos_View_Utils_Def.hpp"

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_MP_VECTOR_CONTIGUOUS_HPP */
