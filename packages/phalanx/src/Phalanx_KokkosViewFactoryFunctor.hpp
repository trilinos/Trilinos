// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_VIEW_FACTORY_FUNCTOR_HPP
#define PHALANX_KOKKOS_VIEW_FACTORY_FUNCTOR_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_Print.hpp"
#include <any>
#include <string>
#include <typeinfo>

namespace PHX {

  template<typename EvalT>
  class KokkosViewFactoryFunctor {
    
    std::unordered_map<std::string,std::any>& fields_;
    const PHX::FieldTag& tag_;
    const std::vector<PHX::index_size_type> extended_dimensions_;
    
  public:
    
    KokkosViewFactoryFunctor(std::unordered_map<std::string,std::any>& fields,
			     const PHX::FieldTag& tag,
			     const std::vector<PHX::index_size_type>& extended_dimensions) :
      fields_(fields),
      tag_(tag),
      extended_dimensions_(extended_dimensions)
    {}
    
    template <typename ScalarT>
    void operator()(ScalarT t) const
    {
      if (tag_.dataTypeInfo() == typeid(t)) {
        using LT = PHX::DataLayout::KokkosLayoutType;
        const auto layout = tag_.dataLayout().kokkosLayout();
        if (layout == LT::Default) {
          PHX::KokkosViewFactory<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> factory;
          fields_[tag_.identifier()] = factory.buildView(tag_,extended_dimensions_);
        }
        else if (layout == LT::Left) {
          PHX::KokkosViewFactory<ScalarT,Kokkos::LayoutLeft,PHX::Device> factory;
          fields_[tag_.identifier()] = factory.buildView(tag_,extended_dimensions_);
        }
        else {
          PHX::KokkosViewFactory<ScalarT,Kokkos::LayoutRight,PHX::Device> factory;
          fields_[tag_.identifier()] = factory.buildView(tag_,extended_dimensions_);
        }
      }
    }
    
  };

}

#endif
