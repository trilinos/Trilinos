// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_TypeAssocMap_hpp__
#define __Panzer_TypeAssocMap_hpp__

#include "Phalanx_TypeList.hpp"

namespace panzer {

/** This class lets you associate evaluation types with
  * a particular value.
  */
template <typename TypesVector,typename ValueType>
class TypeAssocMap {
public:
  typedef TypesVector types_vector;

  TypeAssocMap() 
  {
    constexpr std::size_t sz = PHX::mpl::size_v<TypesVector>;
    mapValues_.resize(sz);
  }

  //! Modify routine 
  template <typename T>
  void set(ValueType v) 
  { 
    constexpr std::size_t idx = PHX::mpl::index_v<T,TypesVector>;
    mapValues_[idx] = v; 
  }

  //! Access routine
  template <typename T>
  ValueType get() const 
  { 
    constexpr std::size_t idx = PHX::mpl::index_v<T,TypesVector>;
    return mapValues_[idx]; 
  }

  template <typename BuilderOpT>
  void buildObjects(const BuilderOpT & builder)
  { PHX::mpl::for_each_type<TypesVector>(BuildObjects<BuilderOpT>(mapValues_,builder)); }

public:

  //! This struct helps will build the values stored in the map
  template <typename BuilderOpT>
  struct BuildObjects {
    std::vector<ValueType> & mv_;
    const BuilderOpT & builder_;
    BuildObjects(std::vector<ValueType> & mv,const BuilderOpT& builder) : mv_(mv), builder_(builder) {}

    template <typename T> void operator()() const {
      constexpr std::size_t idx = PHX::mpl::index_v<T,TypesVector>;
      mv_[idx] = builder_.template build<T>(); 
    }
  };

  std::vector<ValueType> mapValues_;
};
 
}

#endif
