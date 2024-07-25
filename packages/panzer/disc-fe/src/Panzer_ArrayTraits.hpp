// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ARRAY_TRAITS_HPP
#define PANZER_ARRAY_TRAITS_HPP

#include "Kokkos_DynRankView.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

  template<typename Scalar, typename Array> struct ArrayTraits
  {
    typedef typename Array::size_type size_type;
  };
  /*
  // Specialization for Intrepid2::FieldContainer
  template<typename Scalar>
  struct ArrayTraits<Scalar,Kokkos::DynRankView<Scalar,PHX::Device> >
  {
    typedef int size_type;

    // template <typename SubType>
    // struct mod_scalar { typedef Intrepid2::FieldContainer<SubType> array_type; };
    
  };
  */
/*
  // Specialization for MDField
  template<typename Scalar>
  struct ArrayTraits<Scalar,PHX::MDField<Scalar> >
  {
    typedef typename PHX::MDField<Scalar>::size_type size_type;

    // template <typename SubType>
    // struct mod_scalar { typedef PHX::MDField<SubType> array_type; };
  };
*/
}

#endif
