// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Utils.hpp
    \brief  Header function for Intrepid2::Util class and other utility functions.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_UTILS_HPP__
#define __INTREPID2_UTILS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_Types.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_Macros.hpp" // provides some preprocessor values used in definitions of INTREPID2_DEPRECATED, etc.
#include "Kokkos_Random.hpp"

#ifdef HAVE_INTREPID2_SACADO
#include "Kokkos_LayoutNatural.hpp"
#endif

namespace Intrepid2 {

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#define INTREPID2_COMPILE_DEVICE_CODE
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
#define INTREPID2_ENABLE_DEVICE
#endif

#if defined(KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION) \
  && defined(KOKKOS_ENABLE_PRAGMA_IVDEP) \
  && !defined(INTREPID2_COMPILE_DEVICE_CODE)
#define INTREPID2_USE_IVDEP
#endif

  //
  // test macros
  //

#define INTREPID2_TEST_FOR_WARNING(test, msg)                           \
  if (test) {                                                           \
    Kokkos::printf("[Intrepid2] Warning in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
  }

#define INTREPID2_TEST_FOR_EXCEPTION(test, x, msg)                      \
  if (test) {                                                           \
    Kokkos::printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
    throw x(msg);                                                       \
  }

  /// KK: device assert is disabled when NDEBUG is defined which behaves differently 
  ///     from host test.
#ifndef INTREPID2_ENABLE_DEVICE
#define INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(test, x, msg)                              \
  if (test) {                                                                               \
    std::cout << "[Intrepid2] Error in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
    std::cout << "            Test that evaluated to true: " << #test << "\n";              \
    std::cout << "            " << msg << " \n";                                            \
    throw x(msg);                                                                           \
  }
#else
#define INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(test, x, msg)          \
  if (test) {                                                           \
    Kokkos::printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }
#endif
#if defined(INTREPID2_ENABLE_DEBUG) || defined(NDEBUG) || 1
#define INTREPID2_TEST_FOR_ABORT(test, msg)                             \
  if (test) {                                                           \
    Kokkos::printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }
#else
#define INTREPID2_TEST_FOR_ABORT(test, msg) ((void)0)      
#endif
  // check the first error only
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!(info) && (test)) {                                              \
    Kokkos::printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
    info = true;                                                        \
  }
#else  
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!(info) && (test)) {                                              \
    Kokkos::printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    Kokkos::printf("            Test that evaluated to true: %s\n", #test);     \
    Kokkos::printf("            %s \n", msg);                                   \
    info = true ;                                                       \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }
#endif

  /**
   \brief scalar type traits
  */ 
  template<typename T> 
  struct ScalarTraits {
    typedef typename T::scalar_type scalar_type;
  };

  // this is built in types to support
  /**
   \brief Built in support for float
  */ 
  template<>
  struct ScalarTraits<float> {
    typedef float scalar_type;
  };
  /**
   \brief Built in support for double
  */ 
  template<>
  struct ScalarTraits<double> {
    typedef double scalar_type;
  };
  /**
   \brief Built in support for int
  */ 
  template<>
  struct ScalarTraits<int> {
    typedef int scalar_type;
  };
  /**
   \brief Built in support for long int
  */ 
  template<>
  struct ScalarTraits<long int> {
    typedef long int scalar_type;
  };
  /**
   \brief Built in support for long long
  */ 
  template<>
  struct ScalarTraits<long long> {
    typedef long long scalar_type;
  };



  /**
   \brief space overload
  */ 
  template<typename ViewSpaceType, typename UserSpaceType>
  struct ExecSpace {
    typedef UserSpaceType ExecSpaceType;
  };

  /**
   \brief space overload
  */ 
  template<typename ViewSpaceType>
  struct ExecSpace<ViewSpaceType,void> {
    typedef ViewSpaceType ExecSpaceType;
  };


  /**
   \brief layout deduction (temporary meta-function)
  */ 
  template <typename ViewType>
  struct DeduceLayout {
    using input_layout = typename ViewType::array_layout;
    using default_layout = typename ViewType::device_type::execution_space::array_layout;
    using result_layout  =
      typename std::conditional<
        std::is_same< input_layout, Kokkos::LayoutStride >::value,
        default_layout,
        input_layout >::type;
  };


  //
  // utilities device comparible
  //

  // this will be gone 
  template<typename IdxType, typename DimType, typename IterType>
  KOKKOS_FORCEINLINE_FUNCTION
  static void 
  unrollIndex(IdxType &i, IdxType &j, 
              const DimType /* dim0 */,
              const DimType dim1,
              const IterType iter) {
    // left index
    //j = iter/dim0;
    //i = iter%dim0;
    
    // right index
    i = iter/dim1;
    j = iter%dim1;
  }
  
  template<typename IdxType, typename DimType, typename IterType>
  KOKKOS_FORCEINLINE_FUNCTION
  static void 
  unrollIndex(IdxType &i, IdxType &j, IdxType &k, 
              const DimType dim0,
              const DimType dim1,
              const DimType dim2,
              const IterType iter) {
    IdxType tmp;
    
    //unrollIndex(tmp, k, dim0*dim1, dim2, iter);
    //unrollIndex(  i, j, dim0,      dim1,  tmp);
    
    unrollIndex( i, tmp, dim0, dim1*dim2, iter);
    unrollIndex( j, k,   dim1,      dim2,  tmp);
  }
  
  /**
   \brief small utility functions
  */ 
  template<typename T>
  class Util {
  public:
    KOKKOS_FORCEINLINE_FUNCTION
    static T min(const T a, const T b) {
      return (a < b ? a : b);
    }

    KOKKOS_FORCEINLINE_FUNCTION
    static T max(const T a, const T b) {
      return (a > b ? a : b);
    }

    KOKKOS_FORCEINLINE_FUNCTION
    static T abs(const T a) {
      return (a > 0 ? a : T(-a));
    }

  };

  template<typename T>
  KOKKOS_FORCEINLINE_FUNCTION
  static T min(const T &a, const T &b) {
    return (a < b ? a : b);
  }

  template<typename T>
  KOKKOS_FORCEINLINE_FUNCTION
  static T max(const T &a, const T &b) {
    return (a > b ? a : b);
  }

  template<typename T>
  KOKKOS_FORCEINLINE_FUNCTION
  static T abs(const T &a) {
    return (a > 0 ? a : T(-a));
  }

  /**
      \brief functions returning the scalar value.
             for pod types, they return the input object itself.
             of other types it returns the member function val() of the type T.
             For Sacado Fad types it returns the scalar value.
     */

  template<typename T>
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr typename
  std::enable_if< !(std::is_standard_layout<T>::value && std::is_trivial<T>::value), typename ScalarTraits<T>::scalar_type >::type
  get_scalar_value(const T& obj) {return obj.val();}

  template<typename T>
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr typename
  std::enable_if< std::is_standard_layout<T>::value && std::is_trivial<T>::value, typename ScalarTraits<T>::scalar_type >::type
  get_scalar_value(const T& obj){return obj;}


  /**
     \brief specialization of functions for pod types, returning the scalar dimension (1 for pod types) of a view
            these functions are specialized in Sacado and they return the scalar dimension
            (i.e. the dimension of the type w.r.t. the type associated to the pointer type).
    */

  template<typename T, typename ...P>
  KOKKOS_INLINE_FUNCTION
  constexpr typename
  std::enable_if< std::is_standard_layout<T>::value && std::is_trivial<T>::value, unsigned >::type
  dimension_scalar(const Kokkos::DynRankView<T, P...> /* view */) {return 1;}

  template<typename T, typename ...P>
  KOKKOS_INLINE_FUNCTION
  constexpr typename
  std::enable_if< std::is_standard_layout<typename Kokkos::View<T, P...>::value_type>::value && std::is_trivial<typename Kokkos::View<T, P...>::value_type>::value, unsigned >::type
  dimension_scalar(const Kokkos::View<T, P...> /*view*/) {return 1;}

  template<typename T, typename ...P>
  KOKKOS_FORCEINLINE_FUNCTION
  static ordinal_type get_dimension_scalar(const Kokkos::DynRankView<T, P...> &view) {
    return dimension_scalar(view);
  }

  template<typename T, typename ...P>
  KOKKOS_FORCEINLINE_FUNCTION
  static ordinal_type get_dimension_scalar(const Kokkos::View<T, P...> &view) {
    return dimension_scalar(view);
  }
  
  //! \brief Creates and returns a view that matches the provided view in Kokkos Layout.
  //! \param [in] view  - the view to match
  //! \param [in] label - a string label for the view to be created
  //! \param [in] dims  - dimensions to use for the view (the logical dimensions; this method handles adding the derivative dimension required for Fad types).
  //!
  //! This method is particularly useful because we use LayoutStride as the Kokkos Layout in a number of places, and LayoutStride
  //! views cannot be instantiated without also providing stride information.
  //!
  template<class ViewType, class ... DimArgs>
  inline
  Kokkos::DynRankView<typename ViewType::value_type, typename DeduceLayout< ViewType >::result_layout, typename ViewType::device_type >
  getMatchingViewWithLabel(const ViewType &view, const std::string &label, DimArgs... dims)
  {
    using ValueType          = typename ViewType::value_type;
    using ResultLayout       = typename DeduceLayout< ViewType >::result_layout;
    using DeviceType         = typename ViewType::device_type;
    using ViewTypeWithLayout = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;
    
    const bool allocateFadStorage = !(std::is_standard_layout<ValueType>::value && std::is_trivial<ValueType>::value);
    if (!allocateFadStorage)
    {
      return ViewTypeWithLayout(label,dims...);
    }
    else
    {
      const int derivative_dimension = get_dimension_scalar(view);
      return ViewTypeWithLayout(label,dims...,derivative_dimension);
    }
  }

  using std::enable_if_t;

  /**
    \brief Tests whether a class has a member rank.  Used in getFixedRank() method below, which in turn is used in the supports_rank_n helpers.
  */
  template <typename T, typename = void>
  struct has_rank_member : std::false_type{};

  /**
    \brief Tests whether a class has a member rank.  Used in getFixedRank() method below, which in turn is used in the supports_rank_n helpers.
  */
  template <typename T>
  struct has_rank_member<T, decltype((void)T::rank, void())> : std::true_type {};

  static_assert(! has_rank_member<Kokkos::DynRankView<double> >::value, "DynRankView does not have a member rank, so this assert should pass -- if not, something may be wrong with has_rank_member.");
#if KOKKOS_VERSION < 40099
  static_assert(  has_rank_member<Kokkos::View<double*> >::value,        "View has a member rank -- if this assert fails, something may be wrong with has_rank_member.");
#endif

  /**
    \brief \return functor.rank if the functor has a static rank member; returns specified default_value otherwise.
  */
  template<class Functor, ordinal_type default_value>
  constexpr
  enable_if_t<has_rank_member<Functor>::value, ordinal_type>
  getFixedRank()
  {
    return Functor::rank;
  }

  /**
    \brief \return functor.rank if the functor has a static rank member; returns specified default_value otherwise.
  */
  template<class Functor, ordinal_type default_value>
  constexpr
  enable_if_t<!has_rank_member<Functor>::value, ordinal_type>
  getFixedRank()
  {
    return default_value;
  }
 
  /**
    \brief SFINAE helper to detect whether a type supports a 1-integral-argument operator().
  */
  template <typename T>
  class supports_rank_1
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0))>::type );
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,1>() == 1)  };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 2-integral-argument operator().
  */
  template <typename T>
  class supports_rank_2
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,2>() == 2)  };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 3-integral-argument operator().
  */
  template <typename T>
  class supports_rank_3
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = (sizeof(test<T>(0)) == sizeof(char)) && (getFixedRank<T,3>() == 3) };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 4-integral-argument operator().
  */
  template <typename T>
  class supports_rank_4
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0,0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,4>() == 4)  };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 5-integral-argument operator().
  */
  template <typename T>
  class supports_rank_5
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0,0,0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,5>() == 5) };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 6-integral-argument operator().
  */
  template <typename T>
  class supports_rank_6
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0,0,0,0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,6>() == 6)  };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 7-integral-argument operator().
  */
  template <typename T>
  class supports_rank_7
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( typename std::remove_reference<decltype( std::declval<C>().operator()(0,0,0,0,0,0,0))>::type ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) && (getFixedRank<T,7>() == 7) };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a rank-integral-argument operator().
  */
  template <typename T, int rank>
  class supports_rank
  {
  public:
      enum { value = false };
  };

  /**
    \brief SFINAE helper to detect whether a type supports a 1-integral-argument operator().
  */
  template <typename T>
  class supports_rank<T,1>
  {
  public:
      enum { value = supports_rank_1<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 2-integral-argument operator().
  template <typename T>
  class supports_rank<T,2>
  {
  public:
      enum { value = supports_rank_2<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 3-integral-argument operator().
  template <typename T>
  class supports_rank<T,3>
  {
  public:
      enum { value = supports_rank_3<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 4-integral-argument operator().
  template <typename T>
  class supports_rank<T,4>
  {
  public:
      enum { value = supports_rank_4<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 5-integral-argument operator().
  template <typename T>
  class supports_rank<T,5>
  {
  public:
      enum { value = supports_rank_5<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 6-integral-argument operator().
  template <typename T>
  class supports_rank<T,6>
  {
  public:
      enum { value = supports_rank_6<T>::value };
  };

//! SFINAE helper to detect whether a type supports a 7-integral-argument operator().
  template <typename T>
  class supports_rank<T,7>
  {
  public:
      enum { value = supports_rank_7<T>::value };
  };



  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar, int rank>
  struct RankExpander {
    
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,0>
  {
    using value_type = Scalar;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,1>
  {
    using value_type = Scalar*;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,2>
  {
    using value_type = Scalar**;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,3>
  {
    using value_type = Scalar***;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,4>
  {
    using value_type = Scalar****;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,5>
  {
    using value_type = Scalar*****;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,6>
  {
    using value_type = Scalar******;
  };

  /**
    \brief Helper to get Scalar[*+] where the number of *'s matches the given rank.
  */
  template<typename Scalar>
  struct RankExpander<Scalar,7>
  {
    using value_type = Scalar*******;
  };

  // positive checks of supports_rank for Kokkos::DynRankView:
  static_assert(supports_rank<Kokkos::DynRankView<double>, 1>::value, "rank 1 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 2>::value, "rank 2 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 3>::value, "rank 3 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 4>::value, "rank 4 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 5>::value, "rank 5 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 6>::value, "rank 6 check of supports_rank for DynRankView");
  static_assert(supports_rank<Kokkos::DynRankView<double>, 7>::value, "rank 7 check of supports_rank for DynRankView");

  // positive checks of supports_rank for Kokkos::View:
  static_assert(supports_rank<Kokkos::View<double*>,       1>::value, "rank 1 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double**>,      2>::value, "rank 2 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double***>,     3>::value, "rank 3 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double****>,    4>::value, "rank 4 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double*****>,   5>::value, "rank 5 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double******>,  6>::value, "rank 6 check of supports_rank");
  static_assert(supports_rank<Kokkos::View<double*******>, 7>::value, "rank 7 check of supports_rank");

  // negative checks of supports_rank for Kokkos::View:
  static_assert(!supports_rank<Kokkos::View<double*>,       2>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*>,       3>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*>,       4>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*>,       5>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*>,       6>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*>,       7>::value, "rank 1 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      1>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      3>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      4>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      5>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      6>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double**>,      7>::value, "rank 2 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     1>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     2>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     4>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     5>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     6>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double***>,     7>::value, "rank 3 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    1>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    2>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    3>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    5>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    6>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double****>,    7>::value, "rank 4 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   1>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   2>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   3>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   4>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   6>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*****>,   7>::value, "rank 5 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  1>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  2>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  3>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  4>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  5>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double******>,  7>::value, "rank 6 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 1>::value, "rank 7 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 2>::value, "rank 7 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 3>::value, "rank 7 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 4>::value, "rank 7 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 5>::value, "rank 7 check of supports_rank");
  static_assert(!supports_rank<Kokkos::View<double*******>, 6>::value, "rank 7 check of supports_rank");

  /**
    \brief Tests whether a class implements rank().  Used in getFunctorRank() method below; allows us to do one thing for View and another for DynRankView and our custom Functor types.
  */
  template <typename T>
  class has_rank_method
  {
      typedef char one;
      struct two { char x[2]; };

      template <typename C> static one test( decltype( std::declval<C>().rank()  ) ) ;
      template <typename C> static two test(...);

  public:
      enum { value = sizeof(test<T>(0)) == sizeof(char) };
  };

  static_assert(  has_rank_method<Kokkos::DynRankView<double> >::value, "DynRankView implements rank(), so this assert should pass -- if not, something may be wrong with has_rank_method.");
#if KOKKOS_VERSION < 40099
  static_assert(  has_rank_member<Kokkos::View<double*> >::value,        "View has a member rank -- if this assert fails, something may be wrong with has_rank_member.");
#endif

  /**
    \brief \return functor.rank() if the functor implements rank(); functor.rank otherwise.
  */
  template<class Functor>
  enable_if_t<has_rank_method<Functor>::value, unsigned>
  KOKKOS_INLINE_FUNCTION
  getFunctorRank(const Functor &functor)
  {
    return functor.rank();
  }

  /**
    \brief \return functor.rank() if the functor implements rank(); functor.rank otherwise.
  */
  template<class Functor>
  enable_if_t<!has_rank_method<Functor>::value, unsigned>
  KOKKOS_INLINE_FUNCTION
  getFunctorRank(const Functor &functor)
  {
    return functor.rank;
  }

  /**
   \brief Define layout that will allow us to wrap Sacado Scalar objects in Views without copying
   */
#ifdef HAVE_INTREPID2_SACADO
  template <typename ValueType>
  struct NaturalLayoutForType {
    using layout  =
    typename std::conditional<(std::is_standard_layout<ValueType>::value && std::is_trivial<ValueType>::value),
      Kokkos::LayoutLeft, // for POD types, use LayoutLeft
      Kokkos::LayoutNatural<Kokkos::LayoutLeft> >::type; // For FAD types, use LayoutNatural
  };
#else
  template <typename ValueType>
  struct NaturalLayoutForType {
    using layout  = Kokkos::LayoutLeft;
  };
#endif
  
  // define vector sizes for hierarchical parallelism
  const int VECTOR_SIZE = 1;
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD) && defined(INTREPID2_ENABLE_DEVICE)
  const int FAD_VECTOR_SIZE = 32;
#else
  const int FAD_VECTOR_SIZE = 1;
#endif
  
  /**
   \brief Returns a vector size to be used for the provided Scalar type in the context of hierarchically-parallel Kokkos functor execution.
   */
  template<typename Scalar>
  constexpr int getVectorSizeForHierarchicalParallelism()
  {
    return (std::is_standard_layout<Scalar>::value && std::is_trivial<Scalar>::value) ? VECTOR_SIZE : FAD_VECTOR_SIZE;
  }
  
  /**
   \brief Returns the size of the Scalar dimension for the View.  This is 0 for non-AD types.
          This method is useful for sizing scratch storage in hierarchically parallel kernels.
          Whereas get_dimension_scalar() returns 1 for POD types, this returns 0 for POD types.
   */
  template<typename ViewType>
  KOKKOS_INLINE_FUNCTION
  constexpr unsigned getScalarDimensionForView(const ViewType &view)
  {
    return (std::is_standard_layout<typename ViewType::value_type>::value && std::is_trivial<typename ViewType::value_type>::value) ? 0 : get_dimension_scalar(view);
  }
} // end namespace Intrepid2

#endif
