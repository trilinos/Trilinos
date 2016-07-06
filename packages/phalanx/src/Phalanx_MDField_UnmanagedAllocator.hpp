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
    PHX::any memory = PHX::KokkosViewFactory<ScalarT,PHX::Device>::buildView(field.fieldTag(),extra_dims);
    field.setFieldData(memory);
    return field;
  }

} // namespace PHX

#endif
