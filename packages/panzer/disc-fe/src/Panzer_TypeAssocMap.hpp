// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
