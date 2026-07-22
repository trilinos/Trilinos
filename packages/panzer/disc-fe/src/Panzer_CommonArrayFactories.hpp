// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_COMMON_ARRAY_FACTORIES_HPP
#define PANZER_COMMON_ARRAY_FACTORIES_HPP

#include "Kokkos_DynRankView.hpp"
#include "Phalanx_MDField.hpp"

#include <string>

//
// This file contains several common array factories
// useful for build arrays through the BasisValues and
// IntegrationValues classes. In particular these objects
// are used in the <code>setupArrays</code> functions.
// Because these class are used as a template argument the
// types and names used are very specific to the BasisValues
// interface.
//

namespace panzer {
  
  /** Implementation for intrepid field container factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects. Notice in this case the string
    * argument is not used.
    */
  class Intrepid2FieldContainerFactory {
  public:
     template <typename Scalar,typename T0>
     Kokkos::DynRankView<Scalar,PHX::Device> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     Kokkos::DynRankView<Scalar,PHX::Device> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     Kokkos::DynRankView<Scalar,PHX::Device> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     Kokkos::DynRankView<Scalar,PHX::Device> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     Kokkos::DynRankView<Scalar,PHX::Device> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;
  };

  /** Implementation for MDField array factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects.
    */
  class MDFieldArrayFactory {
  public:
     /** Build fields with no prefix, will simply use the string
       * passed into <code>buildArray</code> to name the fields.
       */
     MDFieldArrayFactory() 
      : prefix_(""), allocArray_(false), ddims_(1,0) {}

     /** Build fields with a prefix, will use the string
       * passed into <code>buildArray</code> prefixed with the
       * argument to this constructor to name the fields.
       */
     MDFieldArrayFactory(const std::string & prefix,bool allocArray=false) 
       : prefix_(prefix),allocArray_(allocArray), ddims_(1,0) {}

     /** Build fields with a prefix, will use the string
       * passed into <code>buildArray</code> prefixed with the
       * argument to this constructor to name the fields.
       */
     MDFieldArrayFactory(const std::string & prefix,
                         const std::vector<PHX::index_size_type> & ddims,
                         bool allocArray=false) 
       : prefix_(prefix),allocArray_(allocArray), ddims_(ddims) {}

 
     template <typename Scalar,typename T0>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;

     template <typename Scalar,typename T0>
     PHX::MDField<Scalar,T0> buildStaticArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     PHX::MDField<Scalar,T0,T1> buildStaticArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     PHX::MDField<Scalar,T0,T1,T2> buildStaticArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     PHX::MDField<Scalar,T0,T1,T2,T3> buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     PHX::MDField<Scalar,T0,T1,T2,T3,T4> buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;

  private:
     std::string prefix_;     
     bool allocArray_;
     std::vector<PHX::index_size_type> ddims_;
  };

} // namespace panzer

#include "Panzer_CommonArrayFactories_impl.hpp"

#endif
