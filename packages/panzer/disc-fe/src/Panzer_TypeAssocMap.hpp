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

#include "Sacado_mpl_for_each.hpp"

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
    const int sz = Sacado::mpl::size<TypesVector>::value;
    mapValues_.resize(sz);
  }

  //! Modify routine 
  template <typename T>
  void set(ValueType v) 
  { 
    const int idx = Sacado::mpl::find<TypesVector,T>::value;
    mapValues_[idx] = v; 
  }

  //! Access routine
  template <typename T>
  ValueType get() const 
  { 
    const int idx = Sacado::mpl::find<TypesVector,T>::value;
    return mapValues_[idx]; 
  }

  template <typename BuilderOpT>
  void buildObjects(const BuilderOpT & builder)
  { Sacado::mpl::for_each_no_kokkos<TypesVector>(BuildObjects<BuilderOpT>(mapValues_,builder)); }

public:

  //! This struct helps will build the values stored in the map
  template <typename BuilderOpT>
  struct BuildObjects {
    std::vector<ValueType> & mv_;
    const BuilderOpT & builder_;
    BuildObjects(std::vector<ValueType> & mv,const BuilderOpT& builder) : mv_(mv), builder_(builder) {}

    template <typename T> void operator()(T) const { 
      const int idx = Sacado::mpl::find<TypesVector,T>::value;
      mv_[idx] = builder_.template build<T>(); 
    }
  };

  std::vector<ValueType> mapValues_;
};
 
}

#endif
