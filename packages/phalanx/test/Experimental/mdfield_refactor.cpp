// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <typeinfo>
#include <type_traits>

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// ****************************
// Layouts
// ****************************
template<typename L> struct is_layout : std::false_type {};
template<> struct is_layout<Kokkos::LayoutRight> : std::true_type {};
template<> struct is_layout<Kokkos::LayoutLeft> : std::true_type {};

// ****************************
// Extents
// ****************************
struct C {}; // Cell
struct P {}; // Basis Point
struct Q {}; // Quadrature point
struct D {}; // Dimension

template<typename I> struct is_extent : std::false_type {};
template<> struct is_extent<C> : std::true_type {};
template<> struct is_extent<P> : std::true_type {};
template<> struct is_extent<Q> : std::true_type {};
template<> struct is_extent<D> : std::true_type {};

// ****************************
// Devices
// ****************************
template<typename T> struct is_device : std::false_type {};
template<> struct is_device<PHX::Device> : std::true_type {};

// ****************************
// Rank count
// ****************************
template<typename...Props> struct RankCount;

template<typename...Props> struct RankCount : std::integral_constant<int,RankCount<void,Props...>::value> {};

template<> struct RankCount<void> : std::integral_constant<int,0>{};

template<typename Extent, typename...Props>
struct RankCount<typename std::enable_if<is_extent<Extent>::value>::type,Extent,Props...>
  : std::integral_constant<int,1+RankCount<void,Props...>::value> {};

template<typename NonExtent, typename...Props>
struct RankCount<typename std::enable_if<!is_extent<NonExtent>::value>::type,NonExtent,Props...>
  : std::integral_constant<int,RankCount<void,Props...>::value> {};

// ****************************
// Add pointer (used to construct the static data type)
// ****************************
template<typename Data,int N> struct add_pointer;

template<typename Data,int N> struct add_pointer
{ using type = typename add_pointer<Data*,N-1>::type; };

template<typename Data> struct add_pointer<Data,0>
{ using type = Data; };

// ****************************
// ArrayType
// ****************************
template<typename Scalar,int Rank,typename...Props> struct ArrayType;

// Static rank default
template<typename Scalar,int Rank,typename...Props> struct ArrayType
{
  static_assert(Rank > 0,"Error: Rank of static MDField must be greater than zero!");
  static_assert(Rank < 8,"Error: Rank of static MDField must be less than eight!");
  using data_type = typename add_pointer<Scalar,Rank>::type;
  using array_type = Kokkos::View<data_type,Props...>;
};

// Dynamic rank specialization
template<typename Scalar,typename...Props> struct ArrayType<Scalar,0,Props...>
{
  using data_type = Scalar;
  using array_type = Kokkos::DynRankView<Scalar,Props...>;
};

// ****************************
// FieldTraits
// ****************************
template<typename Scalar, typename... Props> struct FieldTraits;

template<>
struct FieldTraits<void>
{
  using layout = void;
  using device = void;
};

template<typename Extent, typename... Props>
struct FieldTraits< typename std::enable_if<is_extent<Extent>::value>::type, Extent, Props...>
{
  using layout = typename FieldTraits<void,Props...>::layout;
  using device = typename FieldTraits<void,Props...>::device;
};

template<typename Layout, typename... Props>
struct FieldTraits< typename std::enable_if<is_layout<Layout>::value>::type, Layout, Props...>
{
  using layout = Layout;
  using device = typename FieldTraits<void,Props...>::device;
};

template<typename Device, typename... Props>
struct FieldTraits< typename std::enable_if<is_device<Device>::value>::type, Device, Props...>
{
  using layout = void;
  using device = Device;
};

template<typename Scalar, typename... Props>
struct FieldTraits {
  using prop = FieldTraits<void,Props...>;
  static constexpr int rank = RankCount<Props...>::value;
  // This sets defaults if not specified
  using layout = typename std::conditional< !std::is_same<typename prop::layout, void>::value,typename prop::layout, typename PHX::DevLayout<Scalar>::type>::type;
  using device = typename std::conditional< !std::is_same<typename prop::device, void>::value,typename prop::device, PHX::Device>::type;
  using data_type = typename ArrayType<Scalar,rank,layout,device>::data_type;
  using array_type = typename ArrayType<Scalar,rank,layout,device>::array_type;
};

// ****************************
// MDField
// ****************************

template<typename Scalar,typename... Props>
class MDField
{
  // Trick to allow for member method partial specialization on View type (static or dynamic)
  template<int R> struct ViewSpecialization{};

public:
  using traits = FieldTraits<Scalar,Props...>;
  using layout = typename traits::layout;
  using device = typename traits::device;
  using data_type = typename traits::data_type;
  using array_type = typename traits::array_type;
  using size_type = typename traits::array_type::size_type;

  template<typename...Extents>
  MDField(const std::string name,Extents... e)
    : view(name,e...)
  {}

  constexpr bool is_static() const {return (traits::rank != 0);}
  constexpr bool is_dynamic() const {return (traits::rank == 0);}
  constexpr size_type rank() const {return rank(ViewSpecialization<traits::rank>());}


  size_t size() const {return view.size();}
private:
  template<int R> constexpr size_type rank(ViewSpecialization<R>) const {return traits::rank;}
  constexpr size_type rank(ViewSpecialization<0>) const {return view.rank();}


  array_type view;
};

// ****************************
// testing simplifications
// ****************************
template<typename Scalar>
using DefaultLayout = typename PHX::DevLayout<Scalar>::type;

using DefaultDevice = PHX::Device;

// ****************************
// ****************************
// ****************************

TEUCHOS_UNIT_TEST(exp_mdfield_refactor,basic)
{
  // Layouts
  static_assert(is_layout<Kokkos::LayoutLeft>::value,"LayoutLeft broken!");
  static_assert(is_layout<Kokkos::LayoutRight>::value,"LayoutRight broken!");

  // Extents
  static_assert(is_extent<C>::value,"C extent broken!");
  static_assert(is_extent<P>::value,"P extent broken!");
  static_assert(is_extent<Q>::value,"Q extent broken!");
  static_assert(is_extent<D>::value,"D extent broken!");

  // Device
  static_assert(is_device<PHX::Device>::value,"Device broken!");

  // RankCount
  static_assert(RankCount<C,P,D>::value == 3,"RankCount is broken!");
  static_assert(RankCount<C,P,D,PHX::Device>::value == 3,"RankCount is broken!");
  static_assert(RankCount<C,PHX::Device,P,D>::value == 3,"RankCount is broken!");

  // Test add_pointer
  static_assert(std::is_same<double*,typename add_pointer<double,1>::type>::value,"add_pointer is broken");
  static_assert(std::is_same<double**,typename add_pointer<double,2>::type>::value,"add_pointer is broken");
  static_assert(std::is_same<double***,typename add_pointer<double,3>::type>::value,"add_pointer is broken");

  // All defaults
  {
    using ft = FieldTraits<double,C,P>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,DefaultLayout<double>>::value,"default layout is broken!");
    static_assert(std::is_same<typename ft::device,DefaultDevice>::value,"default device is broken!");
    using gold_view = Kokkos::View<double**,DefaultLayout<double>,DefaultDevice>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  // explicit layout, default device
  {
    using ft = FieldTraits<double,C,P,Kokkos::LayoutRight>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,Kokkos::LayoutRight>::value,"explicit layout is broken!");
    static_assert(std::is_same<typename ft::device,DefaultDevice>::value,"default device is broken!");
  }

  // default layout, explicit device
  {
    using ft = FieldTraits<double,C,P,PHX::Device>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,DefaultLayout<double>>::value,"default layout is broken!");
    static_assert(std::is_same<typename ft::device,PHX::Device>::value,"explicit device is broken!");
  }

  // explicit layout, explicit device
  {
    using ft = FieldTraits<double,C,P,Kokkos::LayoutRight,PHX::Device>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,Kokkos::LayoutRight>::value,"explicit layout is broken!");
    static_assert(std::is_same<typename ft::device,PHX::Device>::value,"explicit device is broken!");
    using gold_view = Kokkos::View<double**,Kokkos::LayoutRight,PHX::Device>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  MDField<double,C,P,D> a("a",10,8,3);
  TEST_COMPARE(a.rank(),==,3);
  TEST_COMPARE(a.size(), ==, 10*8*3);
  TEST_COMPARE(a.is_static(), ==, true);
  TEST_COMPARE(a.is_dynamic(), ==, false);

  MDField<double> b("b",10,8,3);
  TEST_COMPARE(b.rank(), ==, 3);
  TEST_COMPARE(b.size(), ==, 10*8*3);
  TEST_COMPARE(b.is_static(), ==, false);
  TEST_COMPARE(b.is_dynamic(), ==, true);
}
