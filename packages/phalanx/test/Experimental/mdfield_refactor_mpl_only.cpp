// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This is a stand-alone toy code for experimenting with template
// metaprogramming wrt MDFields. It is meant to be hand complied
// WITHOUT any Trilinos dependencies using this line:
// g++ --std=c++11 mdfield_refactor_mpl_only.cpp

#include <iostream>
#include <typeinfo>
#include <type_traits>
#include <tuple>

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
// Layouts
// ****************************
struct Left {};
struct Right {};

template<typename L> struct is_layout : std::false_type {};
template<> struct is_layout<Left> : std::true_type {};
template<> struct is_layout<Right> : std::true_type {};

// ****************************
// Devices
// ****************************
struct Serial {};
struct OpenMP {};
struct Cuda {};
template<typename T> struct is_device : std::false_type {};
template<> struct is_device<Serial> : std::true_type {};
template<> struct is_device<OpenMP> : std::true_type {};
template<> struct is_device<Cuda> : std::true_type {};

// ****************************
// Rank count
// ****************************
template<typename...Props> struct RankCount;

template<typename...Props> struct RankCount : std::integral_constant<int,RankCount<void,Props...>::value> {};

template<>
struct RankCount<void>
  : std::integral_constant<int,0>
{};

template<typename Extent, typename...Props>
struct RankCount<typename std::enable_if<is_extent<Extent>::value>::type,Extent,Props...>
  : std::integral_constant<int,1+RankCount<void,Props...>::value>
{};

template<typename NonExtent, typename...Props>
struct RankCount<typename std::enable_if<!is_extent<NonExtent>::value>::type,NonExtent,Props...>
  : std::integral_constant<int,RankCount<void,Props...>::value>
{};

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
  using array_type = std::tuple<data_type,Props...>; // fake for Kokkos::View with static rank
};

// Dynamic rank specialization
template<typename Scalar,typename...Props> struct ArrayType<Scalar,0,Props...>
{
  using data_type = Scalar;
  using array_type = std::tuple<data_type,Props...>; // fake for Kokkos::DynRankView with dynamic rank
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
  // This sets defaults if not specified in template pack
  using layout = typename std::conditional< !std::is_same<typename prop::layout, void>::value,typename prop::layout, Right>::type;
  using device = typename std::conditional< !std::is_same<typename prop::device, void>::value,typename prop::device, Cuda>::type;
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
  using size_type = std::size_t;

  template<typename...Extents>
  MDField(const std::string ,Extents... ) {}

  constexpr bool is_static() const {return (traits::rank != 0);}
  constexpr bool is_dynamic() const {return (traits::rank == 0);}
  constexpr size_type rank() const {return rank(ViewSpecialization<traits::rank>());}

private:
  template<int R> constexpr size_type rank(ViewSpecialization<R>) const {return traits::rank;}
  constexpr size_type rank(ViewSpecialization<0>) const {return 1001;}

  array_type view;
};

// ****************************
// testing simplifications
// ****************************
using DefaultLayout = Right;
using DefaultDevice = Cuda;

// ****************************
// ****************************
// ****************************

int main () {

  // Layouts
  static_assert(is_layout<Left>::value,"LayoutLeft broken!");
  static_assert(is_layout<Right>::value,"LayoutRight broken!");

  // Extents
  static_assert(is_extent<C>::value,"C extent broken!");
  static_assert(is_extent<P>::value,"P extent broken!");
  static_assert(is_extent<Q>::value,"Q extent broken!");
  static_assert(is_extent<D>::value,"D extent broken!");

  // Device
  static_assert(is_device<Serial>::value,"Device broken!");
  static_assert(is_device<OpenMP>::value,"Device broken!");
  static_assert(is_device<Cuda>::value,"Device broken!");

  // RankCount
  static_assert(RankCount<C,P,D>::value == 3,"RankCount is broken!");
  static_assert(RankCount<C,P,D,OpenMP>::value == 3,"RankCount is broken!");
  static_assert(RankCount<C,Cuda,P,D>::value == 3,"RankCount is broken!");
  static_assert(RankCount<C,Cuda,P,Right,D>::value == 3,"RankCount is broken!");

  // Test add_pointer
  static_assert(std::is_same<double*,typename add_pointer<double,1>::type>::value,"add_pointer is broken");
  static_assert(std::is_same<double**,typename add_pointer<double,2>::type>::value,"add_pointer is broken");
  static_assert(std::is_same<double***,typename add_pointer<double,3>::type>::value,"add_pointer is broken");

  // All defaults
  {
    using ft = FieldTraits<double,C,P>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,DefaultLayout>::value,"default layout is broken!");
    static_assert(std::is_same<typename ft::device,DefaultDevice>::value,"default device is broken!");
    using gold_view = std::tuple<double**,DefaultLayout,DefaultDevice>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  // explicit layout, default device
  {
    using ft = FieldTraits<double,C,P,Left>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,Left>::value,"explicit layout is broken!");
    static_assert(std::is_same<typename ft::device,DefaultDevice>::value,"default device is broken!");
    using gold_view = std::tuple<double**,Left,DefaultDevice>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  // default layout, explicit device
  {
    using ft = FieldTraits<double,C,P,Serial>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,DefaultLayout>::value,"default layout is broken!");
    static_assert(std::is_same<typename ft::device,Serial>::value,"explicit device is broken!");
    using gold_view = std::tuple<double**,DefaultLayout,Serial>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  // explicit layout, explicit device
  {
    using ft = FieldTraits<double,C,P,Left,OpenMP>;
    static_assert(ft::rank == 2,"rank is broken!");
    static_assert(std::is_same<typename ft::layout,Left>::value,"explicit layout is broken!");
    static_assert(std::is_same<typename ft::device,OpenMP>::value,"explicit device is broken!");
    using gold_view = std::tuple<double**,Left,OpenMP>;
    static_assert(std::is_same<gold_view,typename ft::array_type>::value,"ArrayType is broken!");
  }

  try {
    MDField<double,C,P,D> a("a",10,8,3);
    if (a.rank() != 3) throw std::runtime_error("incorrect rank");
    if (!a.is_static()) throw std::runtime_error("static check broke");
    if (a.is_dynamic()) throw std::runtime_error("dynamic check broke");
  
    MDField<double> b("b",10,8,3);
    if (b.rank() != 1001) throw std::runtime_error("incorrect rank");
    if (b.is_static()) throw std::runtime_error("static check broke");
    if (!b.is_dynamic()) throw std::runtime_error("dynamic check broke");
  }
  catch (std::exception& e) {
    std::cout << "ERROR: " << e.what() << std::endl; 
    return -1;
  }
  
  return 0;
}
