#include "gtest/gtest.h"
#include "stk_search/HelperTraits.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/IdentProc.hpp"
#include "stk_search/Box.hpp"

namespace {
using BoxIdentType = stk::search::BoxIdent<stk::search::Box<double>, int>;
using BoxIdentProcType = stk::search::BoxIdentProc<stk::search::Box<double>, stk::search::IdentProc<int, int>>;
using IdentIntersectionType = stk::search::IdentIntersection<int, int>;
using IdentProcIntersectionType = stk::search::IdentProcIntersection<stk::search::IdentProc<int, int>, stk::search::IdentProc<int, int>>;

template <typename T, typename U>
struct Foo
{
};

using FooType = Foo<int, int>;
}

TEST(HelperTraits, RemoveCVRef)
{
  static_assert(std::is_same_v<stk::search::remove_cvref_t<int>, int>);
  static_assert(std::is_same_v<stk::search::remove_cvref_t<int&>, int>);
  static_assert(std::is_same_v<stk::search::remove_cvref_t<const int>, int>);
  static_assert(std::is_same_v<stk::search::remove_cvref_t<const int&>, int>);
}

TEST(HelperTraits, ValueTypeOrVoid)
{
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<std::vector<int>>, int>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<std::vector<int>&>, int>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<const std::vector<int>>, int>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<const std::vector<int>&>, int>);

  static_assert(std::is_same_v<stk::search::value_type_or_void_t<int>, void>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<int&>, void>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<const int>, void>);
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<const int&>, void>);
}

TEST(HelperTraits, IsBoxIdent)
{
  static_assert(stk::search::is_box_ident_v<BoxIdentType>);
  static_assert(stk::search::is_box_ident_v<BoxIdentType&>);
  static_assert(stk::search::is_box_ident_v<const BoxIdentType>);
  static_assert(stk::search::is_box_ident_v<const BoxIdentType&>);

  static_assert(!stk::search::is_box_ident_v<int>);
  static_assert(!stk::search::is_box_ident_v<int&>);
  static_assert(!stk::search::is_box_ident_v<const int>);
  static_assert(!stk::search::is_box_ident_v<const int &>);

  static_assert(!stk::search::is_box_ident_v<FooType>);
  static_assert(!stk::search::is_box_ident_v<FooType&>);
  static_assert(!stk::search::is_box_ident_v<const FooType>);
  static_assert(!stk::search::is_box_ident_v<const FooType&>);
}

TEST(HelperTraits, IsBoxIdentProc)
{
  static_assert(stk::search::is_box_ident_proc_v<BoxIdentProcType>);
  static_assert(stk::search::is_box_ident_proc_v<BoxIdentProcType&>);
  static_assert(stk::search::is_box_ident_proc_v<const BoxIdentProcType>);
  static_assert(stk::search::is_box_ident_proc_v<const BoxIdentProcType&>);

  static_assert(!stk::search::is_box_ident_proc_v<int>);
  static_assert(!stk::search::is_box_ident_proc_v<int&>);
  static_assert(!stk::search::is_box_ident_proc_v<const int>);
  static_assert(!stk::search::is_box_ident_proc_v<const int &>);

  static_assert(!stk::search::is_box_ident_proc_v<FooType>);
  static_assert(!stk::search::is_box_ident_proc_v<FooType&>);
  static_assert(!stk::search::is_box_ident_proc_v<const FooType>);
  static_assert(!stk::search::is_box_ident_proc_v<const FooType&>);
}

TEST(HelperTraits, IsBoxIdentContainerSTL)
{
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<std::vector<BoxIdentType>>, BoxIdentType>);
  static_assert(std::is_same_v<stk::search::remove_cvref_t<BoxIdentType>, BoxIdentType>);
  static_assert(stk::search::is_box_ident_container_v<std::vector<BoxIdentType>>);
  static_assert(stk::search::is_box_ident_container_v<std::vector<BoxIdentType>&>);
  static_assert(stk::search::is_box_ident_container_v<const std::vector<BoxIdentType>>);
  static_assert(stk::search::is_box_ident_container_v<const std::vector<BoxIdentType>&>);

  static_assert(!stk::search::is_box_ident_container_v<std::vector<FooType>>);
  static_assert(!stk::search::is_box_ident_container_v<std::vector<FooType>&>);
  static_assert(!stk::search::is_box_ident_container_v<const std::vector<FooType>>);
  static_assert(!stk::search::is_box_ident_container_v<const std::vector<FooType>&>);

  static_assert(!stk::search::is_box_ident_container_v<int>);
  static_assert(!stk::search::is_box_ident_container_v<int&>);
  static_assert(!stk::search::is_box_ident_container_v<const int>);
  static_assert(!stk::search::is_box_ident_container_v<const int&>);
}

TEST(HelperTraits, IsBoxIdentContainerView)
{
  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<BoxIdentType>>);
  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<const BoxIdentType>>);

  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<BoxIdentType>>);
  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<const BoxIdentType>>);

  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<BoxIdentType>&>);
  static_assert(stk::search::is_box_ident_container_v<Kokkos::View<const BoxIdentType>&>);

  static_assert(stk::search::is_box_ident_container_v<const Kokkos::View<BoxIdentType>&>);
  static_assert(stk::search::is_box_ident_container_v<const Kokkos::View<const BoxIdentType>&>);


  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<int>&>);
  static_assert(!stk::search::is_box_ident_container_v<Kokkos::View<const int>&>);

  static_assert(!stk::search::is_box_ident_container_v<const Kokkos::View<int>&>);
  static_assert(!stk::search::is_box_ident_container_v<const Kokkos::View<const int>&>);
}

TEST(HelperTraits, IsBoxIdentProcContainerSTL)
{
  static_assert(std::is_same_v<stk::search::value_type_or_void_t<std::vector<BoxIdentProcType>>, BoxIdentProcType>);
  static_assert(std::is_same_v<stk::search::remove_cvref_t<BoxIdentProcType>, BoxIdentProcType>);
  static_assert(stk::search::is_box_ident_proc_container_v<std::vector<BoxIdentProcType>>);
  static_assert(stk::search::is_box_ident_proc_container_v<std::vector<BoxIdentProcType>&>);
  static_assert(stk::search::is_box_ident_proc_container_v<const std::vector<BoxIdentProcType>>);
  static_assert(stk::search::is_box_ident_proc_container_v<const std::vector<BoxIdentProcType>&>);

  static_assert(!stk::search::is_box_ident_proc_container_v<std::vector<FooType>>);
  static_assert(!stk::search::is_box_ident_proc_container_v<std::vector<FooType>&>);
  static_assert(!stk::search::is_box_ident_proc_container_v<const std::vector<FooType>>);
  static_assert(!stk::search::is_box_ident_proc_container_v<const std::vector<FooType>&>);

  static_assert(!stk::search::is_box_ident_proc_container_v<int>);
  static_assert(!stk::search::is_box_ident_proc_container_v<int&>);
  static_assert(!stk::search::is_box_ident_proc_container_v<const int>);
  static_assert(!stk::search::is_box_ident_proc_container_v<const int&>);
}


TEST(HelperTraits, IsBoxIdentProcContainerView)
{
  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<BoxIdentProcType>>);
  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<const BoxIdentProcType>>);

  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<BoxIdentProcType>>);
  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<const BoxIdentProcType>>);

  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<BoxIdentProcType>&>);
  static_assert(stk::search::is_box_ident_proc_container_v<Kokkos::View<const BoxIdentProcType>&>);

  static_assert(stk::search::is_box_ident_proc_container_v<const Kokkos::View<BoxIdentProcType>&>);
  static_assert(stk::search::is_box_ident_proc_container_v<const Kokkos::View<const BoxIdentProcType>&>);


  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<int>&>);
  static_assert(!stk::search::is_box_ident_proc_container_v<Kokkos::View<const int>&>);

  static_assert(!stk::search::is_box_ident_proc_container_v<const Kokkos::View<int>&>);
  static_assert(!stk::search::is_box_ident_proc_container_v<const Kokkos::View<const int>&>);
}

TEST(HelperTraits, IsIdentIntersection)
{
  static_assert(stk::search::is_ident_intersection_v<IdentIntersectionType>);
  static_assert(stk::search::is_ident_intersection_v<IdentIntersectionType&>);
  static_assert(stk::search::is_ident_intersection_v<const IdentIntersectionType>);
  static_assert(stk::search::is_ident_intersection_v<const IdentIntersectionType&>);

  static_assert(!stk::search::is_ident_intersection_v<FooType>);
  static_assert(!stk::search::is_ident_intersection_v<FooType&>);
  static_assert(!stk::search::is_ident_intersection_v<const FooType>);
  static_assert(!stk::search::is_ident_intersection_v<const FooType&>);
}

TEST(HelperTraits, IsIdentIntersectionContainerSTL)
{
  static_assert(stk::search::is_ident_intersection_container_v<std::vector<IdentIntersectionType>>);
  static_assert(stk::search::is_ident_intersection_container_v<std::vector<IdentIntersectionType>&>);
  static_assert(stk::search::is_ident_intersection_container_v<const std::vector<IdentIntersectionType>>);
  static_assert(stk::search::is_ident_intersection_container_v<const std::vector<IdentIntersectionType>&>);

  static_assert(!stk::search::is_ident_intersection_container_v<std::vector<int>>);
  static_assert(!stk::search::is_ident_intersection_container_v<std::vector<int>&>);
  static_assert(!stk::search::is_ident_intersection_container_v<const std::vector<int>>);
  static_assert(!stk::search::is_ident_intersection_container_v<const std::vector<int>&>);
}

TEST(HelperTraits, IsIdentProcIntersection)
{
  static_assert(stk::search::is_ident_proc_intersection_v<IdentProcIntersectionType>);
  static_assert(stk::search::is_ident_proc_intersection_v<IdentProcIntersectionType&>);
  static_assert(stk::search::is_ident_proc_intersection_v<const IdentProcIntersectionType>);
  static_assert(stk::search::is_ident_proc_intersection_v<const IdentProcIntersectionType&>);

  static_assert(!stk::search::is_ident_proc_intersection_v<FooType>);
  static_assert(!stk::search::is_ident_proc_intersection_v<FooType&>);
  static_assert(!stk::search::is_ident_proc_intersection_v<const FooType>);
  static_assert(!stk::search::is_ident_proc_intersection_v<const FooType&>);
}

TEST(HelperTraits, IsIdentIntersectionContainerView)
{
  static_assert(stk::search::is_ident_intersection_container_v<Kokkos::View<IdentIntersectionType>>);
  static_assert(stk::search::is_ident_intersection_container_v<Kokkos::View<const IdentIntersectionType>>);

  static_assert(stk::search::is_ident_intersection_container_v<Kokkos::View<IdentIntersectionType>&>);
  static_assert(stk::search::is_ident_intersection_container_v<Kokkos::View<const IdentIntersectionType>&>);

  static_assert(stk::search::is_ident_intersection_container_v<const Kokkos::View<IdentIntersectionType>>);
  static_assert(stk::search::is_ident_intersection_container_v<const Kokkos::View<const IdentIntersectionType>>);

  static_assert(stk::search::is_ident_intersection_container_v<const Kokkos::View<IdentIntersectionType>&>);
  static_assert(stk::search::is_ident_intersection_container_v<const Kokkos::View<const IdentIntersectionType>&>);


  static_assert(!stk::search::is_ident_intersection_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_ident_intersection_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_ident_intersection_container_v<Kokkos::View<int>&>);
  static_assert(!stk::search::is_ident_intersection_container_v<Kokkos::View<const int>&>);

  static_assert(!stk::search::is_ident_intersection_container_v<const Kokkos::View<int>>);
  static_assert(!stk::search::is_ident_intersection_container_v<const Kokkos::View<const int>>);

  static_assert(!stk::search::is_ident_intersection_container_v<const Kokkos::View<int>&>);
  static_assert(!stk::search::is_ident_intersection_container_v<const Kokkos::View<const int>&>);
}

TEST(HelperTraits, IsIdentProcIntersectionContainerSTL)
{
  static_assert(stk::search::is_ident_proc_intersection_container_v<std::vector<IdentProcIntersectionType>>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<std::vector<IdentProcIntersectionType>&>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<const std::vector<IdentProcIntersectionType>>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<const std::vector<IdentProcIntersectionType>&>);

  static_assert(!stk::search::is_ident_proc_intersection_container_v<std::vector<int>>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<std::vector<int>&>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<const std::vector<int>>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<const std::vector<int>&>);
}

TEST(HelperTraits, IsIdentProcIntersectionContainerView)
{
  static_assert(stk::search::is_ident_proc_intersection_container_v<Kokkos::View<IdentProcIntersectionType>>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<Kokkos::View<const IdentProcIntersectionType>>);

  static_assert(stk::search::is_ident_proc_intersection_container_v<Kokkos::View<IdentProcIntersectionType>&>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<Kokkos::View<const IdentProcIntersectionType>&>);

  static_assert(stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<IdentProcIntersectionType>>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<const IdentProcIntersectionType>>);

  static_assert(stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<IdentProcIntersectionType>&>);
  static_assert(stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<const IdentProcIntersectionType>&>);


  static_assert(!stk::search::is_ident_proc_intersection_container_v<Kokkos::View<int>>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<Kokkos::View<const int>>);

  static_assert(!stk::search::is_ident_proc_intersection_container_v<Kokkos::View<int>&>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<Kokkos::View<const int>&>);

  static_assert(!stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<int>>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<const int>>);

  static_assert(!stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<int>&>);
  static_assert(!stk::search::is_ident_proc_intersection_container_v<const Kokkos::View<const int>&>);
}

TEST(HelperTraits, IsModifiable)
{
  static_assert(stk::search::is_modifiable_v<int>);
  static_assert(stk::search::is_modifiable_v<int&>);
  static_assert(!stk::search::is_modifiable_v<const int>);
  static_assert(!stk::search::is_modifiable_v<const int&>);
}

TEST(HelperTraits, IsModifiableView)
{
  static_assert(stk::search::is_modifiable_view_v<Kokkos::View<int>>);
  static_assert(stk::search::is_modifiable_view_v<Kokkos::View<int>&>);
  static_assert(stk::search::is_modifiable_view_v<const Kokkos::View<int>>);
  static_assert(stk::search::is_modifiable_view_v<const Kokkos::View<int>&>);

  static_assert(!stk::search::is_modifiable_view_v<Kokkos::View<const int>>);
  static_assert(!stk::search::is_modifiable_view_v<Kokkos::View<const int>&>);
  static_assert(!stk::search::is_modifiable_view_v<const Kokkos::View<const int>>);
  static_assert(!stk::search::is_modifiable_view_v<const Kokkos::View<const int>&>);
}

TEST(HelperTraits, ViewUsableFrom)
{
  stk::search::check_view_is_usable_from<Kokkos::DefaultHostExecutionSpace, Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace>>();
  stk::search::check_view_is_usable_from<Kokkos::DefaultHostExecutionSpace, Kokkos::View<int*, Kokkos::HostSpace>>();
  stk::search::check_view_is_usable_from<Kokkos::DefaultHostExecutionSpace, Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace>>();
  stk::search::check_view_is_usable_from<Kokkos::DefaultHostExecutionSpace, Kokkos::View<int**, Kokkos::LayoutRight, Kokkos::HostSpace>>();

#ifdef KOKKOS_ENABLE_CUDA
  stk::search::check_view_is_usable_from<Kokkos::Cuda, Kokkos::View<int*, Kokkos::Cuda>>();
  stk::search::check_view_is_usable_from<Kokkos::Cuda, Kokkos::View<int*, Kokkos::CudaSpace>>();
  stk::search::check_view_is_usable_from<Kokkos::Cuda, Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::CudaSpace>>();
  stk::search::check_view_is_usable_from<Kokkos::Cuda, Kokkos::View<int**, Kokkos::LayoutRight, Kokkos::CudaSpace>>();
#endif
}
