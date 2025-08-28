// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_PRINT_VALUES_HPP
#define PHALANX_PRINT_VALUES_HPP

#include <ostream>
#include <string>
#include <type_traits>
#include <Kokkos_Core.hpp>
#include <Phalanx_MDField.hpp>
#include <Phalanx_Field.hpp>

namespace PHX {

  // Specializations to check for type of array (Kokkos already supports is_view for View)
  template<typename View> struct is_dyn_rank_view : std::false_type {};
  template<typename T, typename ...Props> struct is_dyn_rank_view<Kokkos::DynRankView<T,Props...>> : std::true_type {};
  template<typename T, typename ...Props> struct is_dyn_rank_view<const Kokkos::DynRankView<T,Props...>> : std::true_type {};

  template<typename View> struct is_mdfield : std::false_type {};
  template<typename T, typename ...Props> struct is_mdfield<PHX::MDField<T,Props...>> : std::true_type {};
  template<typename T, typename ...Props> struct is_mdfield<const PHX::MDField<T,Props...>> : std::true_type {};

  template<typename View> struct is_field : std::false_type {};
  template<typename T, int Rank, typename Layout> struct is_field<PHX::Field<T,Rank,Layout>> : std::true_type {};
  template<typename T, int Rank, typename Layout> struct is_field<const PHX::Field<T,Rank,Layout>> : std::true_type {};

  /// Function to print values of Kokkos::View, Kokkos::DynRankView, PHX::MDField and PHX::Field 
  template<typename ViewType>
  void printValues(const ViewType& view_device, std::ostream& os, std::string prefix = "")
  {
    static_assert(!std::is_pointer<ViewType>::value);
    
    if constexpr (is_dyn_rank_view<ViewType>::value) {
      typename ViewType::host_mirror_type v = Kokkos::create_mirror(view_device);
      Kokkos::deep_copy(v,view_device);
      const auto rank = v.rank();
      if (rank == 1) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          os << prefix << "(" << i << ") = " << v(i) << std::endl;
      }
      else if (rank == 2) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            os << prefix << "(" << i << "," << j << ") = " << v(i,j) << std::endl;
      }
      else if (rank == 3) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              os << prefix << "(" << i << "," << j << "," << k << ") = " << v(i,j,k) << std::endl;
      }
      else if (rank == 4) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                os << prefix << "(" << i << "," << j << "," << k << "," << m << ") = "
                   << v(i,j,k,m) << std::endl;
      }
      else if (rank == 5) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                     << n << ") = "
                     << v(i,j,k,m,n) << std::endl;
      }
      else if (rank == 6) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  for (std::size_t p=0; p < v.extent(5); ++p)
                    os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                       << n << "," << p << ") = "
                       << v(i,j,k,m,n,p) << std::endl;
      }
      else if (rank == 7) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  for (std::size_t p=0; p < v.extent(5); ++p)
                    for (std::size_t q=0; q < v.extent(6); ++q)
                      os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                         << n << "," << p << "," << q << ") = "
                         << v(i,j,k,m,n,p,q) << std::endl;
      }
    }
    else if constexpr (Kokkos::is_view<ViewType>::value){ // Static rank View
      typename ViewType::host_mirror_type v = Kokkos::create_mirror(view_device);
      Kokkos::deep_copy(v,view_device);
      if constexpr (v.rank()==1) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          os << prefix << "(" << i << ") = " << v(i) << std::endl;
      }
      else if constexpr (v.rank()==2) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            os << prefix << "(" << i << "," << j << ") = " << v(i,j) << std::endl;
      }
      else if constexpr (v.rank()==3) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              os << prefix << "(" << i << "," << j << "," << k << ") = " << v(i,j,k) << std::endl;
      }
      else if constexpr (v.rank()==4) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                os << prefix << "(" << i << "," << j << "," << k << "," << m << ") = "
                   << v(i,j,k,m) << std::endl;
      }
      else if constexpr (v.rank()==5) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                     << n << ") = "
                     << v(i,j,k,m,n) << std::endl;
      }
      else if constexpr (v.rank()==6) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  for (std::size_t p=0; p < v.extent(5); ++p)
                    os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                       << n << "," << p << ") = "
                       << v(i,j,k,m,n,p) << std::endl;
      }
      else if constexpr (v.rank()==7) {
        for (std::size_t i=0; i < v.extent(0); ++i)
          for (std::size_t j=0; j < v.extent(1); ++j)
            for (std::size_t k=0; k < v.extent(2); ++k)
              for (std::size_t m=0; m < v.extent(3); ++m)
                for (std::size_t n=0; n < v.extent(4); ++n)
                  for (std::size_t p=0; p < v.extent(5); ++p)
                    for (std::size_t q=0; q < v.extent(6); ++q)
                      os << prefix << "(" << i << "," << j << "," << k << "," << m << ","
                         << n << "," << p << "," << q << ") = "
                         << v(i,j,k,m,n,p,q) << std::endl;
      }
    }
    else if constexpr (is_mdfield<ViewType>::value) {
      if constexpr (ViewType::traits::rank == 0) {
        PHX::printValues(view_device.get_view(),os,prefix);
      }
      else {
        PHX::printValues(view_device.get_static_view(),os,prefix);
      }
    }
    else if constexpr (is_field<ViewType>::value) {
      PHX::printValues(view_device.get_static_view(),os,prefix);
    }
    else {
      static_assert(!Kokkos::is_view<ViewType>::value &&
                    !is_dyn_rank_view<ViewType>::value &&
                    !is_mdfield<ViewType>::value &&
                    !is_field<ViewType>::value,
                    "Error: PHX::printValues() only supports Kokkos::View, Kokkos::DynRankView, PHX::MDField and PHX::Field!");
    }
  }
}

#endif
