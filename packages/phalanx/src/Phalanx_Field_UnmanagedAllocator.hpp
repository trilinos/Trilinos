#ifndef PHALANX_FIELD_UNMANAGED_ALLOCATOR_HPP
#define PHALANX_FIELD_UNMANAGED_ALLOCATOR_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_KokkosViewFactory.hpp"

namespace PHX {

  /*! \brief Allocates an UNMANAGED compiletime Field from a PHX::DataLayout. */
  template<typename ScalarT,int Rank>
  PHX::Field<ScalarT,Rank>
  allocateUnmanagedField(const std::string& name,
                         const Teuchos::RCP<PHX::DataLayout>& layout,
                         const std::vector<PHX::index_size_type>& extra_dims = std::vector<PHX::index_size_type>(0))
  {
    PHX::Field<ScalarT,Rank> field(name,layout);
    PHX::any memory = PHX::KokkosViewFactory<ScalarT,PHX::Device>::buildView(field.fieldTag(),extra_dims);
    field.setFieldData(memory);
    return field;
  }

} // namespace PHX

#endif
