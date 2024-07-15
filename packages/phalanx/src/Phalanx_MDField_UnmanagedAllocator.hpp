// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_MDFIELD_UNMANAGED_ALLOCATOR_HPP
#define PHALANX_MDFIELD_UNMANAGED_ALLOCATOR_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosViewFactory.hpp"

namespace PHX {

  /*! \brief Allocates an UNMANAGED compiletime MDField from a PHX::DataLayout. */
  template<typename ScalarT,typename ... DimensionPack>
  PHX::MDField<ScalarT,DimensionPack...>
  allocateUnmanagedMDField(const std::string& name,
                           const Teuchos::RCP<PHX::DataLayout>& layout,
                           const std::vector<PHX::index_size_type>& extra_dims = std::vector<PHX::index_size_type>(0))
  {
    PHX::MDField<ScalarT,DimensionPack...> field(name,layout);
    std::any memory = PHX::KokkosViewFactory<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>::buildView(field.fieldTag(),extra_dims);
    field.setFieldData(memory);
    return field;
  }

} // namespace PHX

#endif
