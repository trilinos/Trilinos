#ifndef PHX_DATA_CONTAINER_BASE_HPP
#define PHX_DATA_CONTAINER_BASE_HPP

#include <typeinfo>
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"

namespace PHX {

  template<typename Traits>
  class DataContainerBase {

  public:

    DataContainerBase();

    virtual ~DataContainerBase();

    virtual void allocateField(const Teuchos::RCP<PHX::FieldTag>& v,
			       std::size_t max_num_cells,
			       typename Traits::Allocator& a) = 0;
    
    virtual const std::type_info& dataTypeInfo() const = 0; 

    virtual std::size_t getSizeOfDataType() const = 0;

    virtual void print(std::ostream& os) const = 0;
    
  };

  template<typename Traits>
  std::ostream& 
  operator<<(std::ostream& os, 
	     const PHX::DataContainerBase<Traits>& dc);

}

#include "Phalanx_DataContainer_Base_Def.hpp"

#endif 
